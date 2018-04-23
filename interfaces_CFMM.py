#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 11:12:06 2017

@author: akuurstr
"""
from scipy import ndimage
import traits.api as traits
import nibabel as nib
import numpy as np
import os,pyQSM.frequencyEstimate
from nipype.interfaces.base import CommandLineInputSpec,BaseInterface,TraitedSpec,InputMultiPath,File,CommandLine,OutputMultiPath
from nipype.interfaces.traits_extension import isdefined
import pyQSM.calculateReliability as cr
import pyQSM.dipoleInversion_Liu2012
import os.path
import json
from scipy.optimize import curve_fit,leastsq
from nipype.interfaces.utility import Function

matlab_scripts_loc=os.path.dirname(os.path.realpath(__file__))+"/matlab_scripts"
#mcr_location='/usr/local/MATLAB/R2016b/' #local testing
#r2_script_location=os.path.dirname(os.path.realpath(__file__))+'/CalcR2Star.py' #local testing
mcr_location='/opt/mcr/v91'
r2_script_location='/code/CalcR2Star.py'

#==============================================================================
#SiemensPhasePreprocess interface is for converting siemens dicom to radians
#and any making sure the phase is in right handed system
#==============================================================================
def siemens2rad(phase_img_siemens):
        #This function converts an n dimensional array from Siemens phase units to
        #radians.  Siemens employs the following phase convention: 
        #S = 2048[(R/pi) + 1], where 'S' is the phase in Siemens units (i.e. what
        #is returned by the scanner) and 'R' is the phase in radians. The range of
        #any phase image, in radians, will be -pi to +pi.
        #The output, 'phase_img_rad', is of class double.
        
        phase_img_siemens = phase_img_siemens.astype('float');
        phase_img_rad = (phase_img_siemens/2048 - 1)*np.pi;
        return phase_img_rad

class SiemensPhasePreprocessInputSpec(CommandLineInputSpec):
    infiles = InputMultiPath(File(exists=True), desc='A list of phase echoes',
                                 copyFile=False, mandatory=True)
    
class SiemensPhasePreprocessOutputSpec(TraitedSpec):    
    outfiles = OutputMultiPath(File(exists=True), desc='A list of vendor-specific processed phase echoes',
                                 copyFile=False, mandatory=True)
     
class SiemensPhasePreprocess(BaseInterface):
    input_spec = SiemensPhasePreprocessInputSpec
    output_spec = SiemensPhasePreprocessOutputSpec    
    
    def _run_interface(self, runtime):
        infiles=self.inputs.infiles
        self.outfilenames=[]
                
        for f in infiles:
            imgobj=nib.load(f)
            img=imgobj.get_data()
            img=siemens2rad(img)
            filename,ext=os.path.splitext(f)
            ext2=''
            if ext=='.gz':
                filename,ext2=os.path.splitext(filename)                
            newfilename=os.path.abspath(os.path.basename(filename)+'_processed'+ext2+ext)
            niftifile=nib.Nifti1Pair(img,imgobj.affine)
            nib.save(niftifile,newfilename)
            self.outfilenames.append(newfilename)              
        return runtime
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['outfiles'] = self.outfilenames
        return outputs 

"""
#quick and dirty but no control over input/output spec
def SiemensPhasePreprocess_fn(infiles):
    import nibabel as nib
    import numpy as np
    import os
    def siemens2rad(phase_img_siemens):      
        phase_img_siemens = phase_img_siemens.astype('float');
        phase_img_rad = (phase_img_siemens/2048 - 1)*np.pi;
        return phase_img_rad    
    outfilenames=[]
    for f in infiles:
        imgobj=nib.load(f)
        img=imgobj.get_data()
        img=siemens2rad(img)
        filename,ext=os.path.splitext(f)
        newfilename=os.path.abspath(os.path.basename(filename)+'_processed'+ext)            
        niftifile=nib.Nifti1Pair(img,imgobj.affine)
        nib.save(niftifile,newfilename)
        outfilenames.append(newfilename)
    return outfilenames    
SiemensPhasePreprocess = Function(input_names=["infiles"],
                             output_names=["outfiles"],
                             function=SiemensPhasePreprocess_fn)
#"""    
#==============================================================================    

#==============================================================================
#VarianPhasePreprocess makes sure the phase is in the correct handed system
#==============================================================================
class VarianPhasePreprocessInputSpec(CommandLineInputSpec):
    infiles = InputMultiPath(File(exists=True), desc='A list of phase echoes',
                                 copyFile=False, mandatory=True)
    
class VarianPhasePreprocessOutputSpec(TraitedSpec):    
    outfiles = OutputMultiPath(File(exists=True), desc='A list of vendor-specific processed phase echoes',
                                 copyFile=False, mandatory=True)
     
class VarianPhasePreprocess(BaseInterface):
    input_spec = VarianPhasePreprocessInputSpec
    output_spec = VarianPhasePreprocessOutputSpec    
    
    def _run_interface(self, runtime):
        infiles=self.inputs.infiles
        self.outfilenames=[]
                
        for f in infiles:
            imgobj=nib.load(f)
            img=imgobj.get_data()
            img=-img
            filename,ext=os.path.splitext(f)
            newfilename=os.path.abspath(os.path.basename(filename)+'_processed'+ext)
            niftifile=nib.Nifti1Pair(img,imgobj.affine)
            nib.save(niftifile,newfilename)
            self.outfilenames.append(newfilename)              
        return runtime
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['outfiles'] = self.outfilenames
        return outputs 
#==============================================================================    

#==============================================================================
#calcR2star Interface is for calculating r2* from complex multi-echo data.
#it is written as a python interface and uses multiprocessing to parallelize across
#spatial voxels. in python2 each worker is created using fork, and the address space is 
#copied on write. since we only read a pixel's multi-echo data to compute
#r2*, the large input data should not actually be copied for each worker.
#however, compute canada servers count the copy on write requirements towards
#your total memory usage and this method quickly runs out of memory.

#calcR2star_cmd Interface calls a r2* script using python3.  The r2* script is
#similar to the calcR2star Interface code, but parallelizes using python3
#multiprocessing with 'spawn' instead of fork to properly report memory usage.
#==============================================================================
def model_monoexponential_R2Star(TE,A,B,R2star,f): 
    #do not use mag/phase for scaling parameter A!!! use re/imag
    return (A+1j*B)*np.exp(-TE*R2star+1j*f*TE)
    
def residuals(params,TE,data):
    diff=model_monoexponential_R2Star(TE,*params)-data
    diffconcat=np.empty(diff.size*2,dtype='float')
    diffconcat[0:diffconcat.size:2]= diff.real
    diffconcat[1:diffconcat.size:2]= diff.imag   
    return diffconcat

#first implementation, serially processed
#def calcR2starFn(img,mask,TE):     
#    img_shape=img.shape
#    img=img/(img[...,0]*mask).max() #useful if you want to use ftol to exit the curve_fit
#    img=img.reshape(np.prod(img_shape[:-1]),img_shape[-1])    
#    mask=mask.reshape(np.prod(img_shape[:-1]))
#    #Amplitude=np.empty(np.prod(img_shape[:-1]))
#    R2star=np.empty(np.prod(img_shape[:-1]))
#    goodnessOfFit_r2=np.empty(np.prod(img_shape[:-1]))
#    """    
#    r2InitEstimate=-np.log(np.abs(img)[...,1]/np.abs(img)[...,0])/(TE[1]-TE[0])
#    r2InitEstimate=r2InitEstimate*(r2InitEstimate>0)
#    AInitEstimate=np.abs(img)[...,0]/np.exp(-TE[0]*r2InitEstimate)
#    """
#    A=np.array([np.ones_like(TE),-TE]).T   
#    y=np.log(img.reshape(np.prod(img.shape[:-1]),img.shape[-1])+1e-20).T 
#    x=np.dot(np.linalg.pinv(A),y)         
#
#    for indx in range(img.shape[0]):        
#        if mask[indx]:
#            try:                
#                scale,exponent_factor=x[:,indx]                    
#                x0=(scale.real,scale.imag,exponent_factor.real,exponent_factor.imag)                
#                lstsqopt=leastsq(residuals,x0,maxfev=10000,args=(TE,img[indx,:]),full_output=True)
#                ss_err=(lstsqopt[2]['fvec']**2).sum()
#                ss_tot=(np.abs(img[indx]-img[indx].mean())**2).sum()
#                rsquared=1-(ss_err/ss_tot)
#                R2star[indx]=lstsqopt[0][2]
#                goodnessOfFit_r2[indx]=rsquared
#            except:
#                R2star[indx]=np.nan   
#                goodnessOfFit_r2[indx]=np.nan
#        else:
#            R2star[indx]=0                
#    R2star=R2star.reshape(img_shape[:-1])
#    goodnessOfFit_r2=goodnessOfFit_r2.reshape(img_shape[:-1])
#    negMask=R2star<0
#    nanMask=np.isnan(R2star)
#    R2star=R2star*np.invert(negMask)*np.invert(nanMask)
#    return R2star,goodnessOfFit_r2,negMask,nanMask


def calc_R2star_parallelizable_unit(img,TE,initial_guess,indx):
    initial_scale,initial_exponent_factor=initial_guess[:,indx]
    x0=(initial_scale.real,initial_scale.imag,initial_exponent_factor.real,initial_exponent_factor.imag)
    #x0=(scale.real,scale.imag,exponent_factor.real)
    lstsqopt=leastsq(residuals,x0,args=(TE,img[indx]),maxfev=10000,full_output=True)
 
    ss_err=(lstsqopt[2]['fvec']**2).sum()
    ss_tot=(np.abs(img[indx]-img[indx].mean())**2).sum()
    rsquared=1-(ss_err/ss_tot)

    if lstsqopt[-1] in (1,2,3,4):
        return lstsqopt[0][2],rsquared
    else:        
        return np.nan,np.inf

           
def calc_R2star_fn(img,mask,TE):
    img_shape=img.shape
    img=img/(img[...,0]*mask).max() #useful if you want to use ftol to exit the curve_fit
    img=img.reshape(np.prod(img_shape[:-1]),img_shape[-1])    
    
    #initial variable guess
    A=np.array([np.ones_like(TE),-TE]).T   
    y=np.log(img.reshape(np.prod(img.shape[:-1]),img.shape[-1])+1e-20).T 
    x=np.dot(np.linalg.pinv(A),y)    
    
    indices=np.arange(np.prod(img_shape[:-1]))
    indices=indices.reshape(img_shape[:-1])    
    indices=indices[mask.astype('bool')]
            
    from multiprocessing import Pool        
    pool = Pool()
    
    from functools import partial
    calc_R2star_parallelizable_unit_partial = partial(calc_R2star_parallelizable_unit, img,TE,x)    
    result=np.array(pool.map(calc_R2star_parallelizable_unit_partial, indices))
    pool.close()
    pool.join()

    R2star=result[:,0]
    goodness_of_fit_r2=result[:,1]
    
    tmp=np.zeros(np.prod(img_shape[:-1]))
    tmp[indices]=R2star
    tmp=tmp.reshape(img_shape[:-1])
    R2star=tmp
    
    tmp=np.zeros(np.prod(img_shape[:-1]))
    tmp[indices]=goodness_of_fit_r2
    tmp=tmp.reshape(img_shape[:-1])
    goodness_of_fit_r2=tmp
    
    neg_mask=R2star<0
    nan_mask=np.isnan(R2star)
    R2star=R2star*np.invert(neg_mask)*np.invert(nan_mask)
    return R2star,goodness_of_fit_r2,neg_mask,nan_mask
            
class CalcR2StarInputSpec(CommandLineInputSpec):
    mag = InputMultiPath(File(exists=True), desc='A list of magnitude echoes',
                          argstr="--mag_filename_list %s", position=1,
                          copyFile=False, mandatory=True)
    phase = InputMultiPath(File(exists=True), desc='A list of phase echoes',
                           argstr="--phase_filename_list %s", position=2,
                           copyFile=False, mandatory=True)
    json = InputMultiPath(File(exists=True), desc='A list of json files containing the tes',
                      argstr="--json_filename_list %s", position=3,
                      copyFile=False, mandatory=True)        
    freq_loc = File(exists=True, desc='Input freq filename',      
                   argstr="--freq_filename %s", position=4,                           
                   copyFile=False, mandatory=True)    
    mask = File(desc = 'Mask (3D) where to calculate phase reliability',
                argstr="--mask_filename %s", position=5,
                copyFile=False,mandatory=True)
    R2star = File(desc="Output R2star file",mandatory=True,
                  argstr="--r2star_filename %s", position=6,)    
    neg_mask = File(desc="Mask of pixels that fit to a negative R2star and were set to 0",
                   mandatory=True,argstr="--neg_mask_filename %s", position=7,)
    nan_mask = File(desc="Mask of pixels where fitting failed and were set to 0",
                   mandatory=True,argstr="--nan_mask_filename %s", position=8,)
    
class CalcR2StarOutputSpec(TraitedSpec):    
    R2star = File(desc="Output file containing the estimated R2star",
                   exists=True)
    R2star_fit = File(desc="Output file containing the R-squared measure of goodness of fit",
                   exists=True)
    neg_mask = File(desc="Output file containing the mask of pixels that fit to a negative R2star and were set to 0",
                   exists=True)
    nan_mask = File(desc="Output file containing the mask of pixels where fitting failed and were set to 0",
                   exists=True)
    
class CalcR2Star(BaseInterface):
    input_spec = CalcR2StarInputSpec
    output_spec = CalcR2StarOutputSpec    
    
    def _run_interface(self, runtime):
        mag_filename_list=self.inputs.mag
        phase_filename_list=self.inputs.phase
        freq_filename=self.inputs.freq_loc
        mask_filename=self.inputs.mask
        json_filename_list=self.inputs.json        
                
        R2star_filename=self._list_outputs()['R2star']
        neg_mask_filename=self._list_outputs()['neg_mask']        
        nan_mask_filename=self._list_outputs()['nan_mask']        
        
        mag_img_obj=nib.load(mag_filename_list[0])
        shape=mag_img_obj.get_shape()
        
        complex_img=np.empty(shape+(len(mag_filename_list),),dtype='complex')        
        te=np.empty(len(json_filename_list))
        count=0
        
        for magloc,phloc,jsonloc in zip(mag_filename_list,phase_filename_list,json_filename_list):
            complex_img[...,count]=(nib.load(magloc).get_data()).astype('float')*np.exp(1j*(nib.load(phloc).get_data())) 
            with open(jsonloc) as f:
                te[count]=json.load(f)['EchoTime']
            count+=1        
        freq=nib.load(freq_filename).get_data()
        mask=nib.load(mask_filename).get_data()
        #to avoid temporal wraps, we can remove the phase evolution of each echo
        #we retain the noise useful for the complex fit
        data_reduced_phase=complex_img*np.exp(-1j*freq[...,np.newaxis]*te) 
        R2star_img,goodness_of_fit,neg_mask,nan_mask=calc_R2star_fn(data_reduced_phase,mask,te)
        neg_mask=neg_mask.astype('int8')
        nan_mask=nan_mask.astype('int8')
        
        niftifile=nib.Nifti1Pair(R2star_img,mag_img_obj.affine)
        nib.save(niftifile,R2star_filename)
        niftifile=nib.Nifti1Pair(goodness_of_fit,mag_img_obj.affine)
        nib.save(niftifile,os.path.splitext(R2star_filename)[0]+'_fit'+os.path.splitext(R2star_filename)[1])
        niftifile=nib.Nifti1Pair(neg_mask,mag_img_obj.affine)
        nib.save(niftifile,neg_mask_filename)
        niftifile=nib.Nifti1Pair(nan_mask,mag_img_obj.affine)
        nib.save(niftifile,nan_mask_filename)
        
        return runtime
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['R2star'] = os.path.abspath(self.inputs.R2star)
        outputs['R2star_fit'] = os.path.abspath(os.path.splitext(self.inputs.R2star)[0]+'_fit'+os.path.splitext(self.inputs.R2star)[1])
        outputs['neg_mask'] = os.path.abspath(self.inputs.neg_mask)
        outputs['nan_mask'] = os.path.abspath(self.inputs.nan_mask)
        return outputs 
    
class CalcR2Star_cmd(CommandLine):
    input_spec = CalcR2StarInputSpec
    output_spec = CalcR2StarOutputSpec  
    cmd = 'python3 '+r2_script_location
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['R2star'] = os.path.abspath(self.inputs.R2star)
        #extension should start at leading dot instead of last dot
        #use splitext in reverse
        ext,name=os.path.splitext(self.inputs.R2star[::-1])
        name=name[:0:-1]
        ext='.'+ext[::-1]
        outputs['R2star_fit'] = os.path.abspath(name+'_fit'+ext)
        outputs['neg_mask'] = os.path.abspath(self.inputs.neg_mask)
        outputs['nan_mask'] = os.path.abspath(self.inputs.nan_mask)
        return outputs 
#==============================================================================    

#==============================================================================
# ImHistMatch interface uses matlab to histogram match a list of images to a reference
#==============================================================================
    
class ImHistMatchInputSpec(CommandLineInputSpec):

    in_file = File(exists=True, desc='Input filename',
                                 argstr='%s', position=1,
                                 copyFile=False, mandatory=True)
        
    ref = File(exists=True, desc = 'Reference filename', argstr='%s',
                  position=2,copyFile=False,mandatory=True)  
    
    result_filename = File(desc="Output LFS filename",argstr='%s',position=3,               
                   mandatory=True)   
    
    quit_matlab = traits.String(';quit;"',desc='needed to quit matlab',argstr='%s',
                      position=4,mandatory=False,usedefault=True)   
    
class ImHistMatchOutputSpec(TraitedSpec):    
    result_filename = File(desc="Histogram matched image",exists=True)
    
    
class ImHistMatch(CommandLine):
    input_spec = ImHistMatchInputSpec
    output_spec = ImHistMatchOutputSpec
    cmd = 'matlab -nodisplay -nosplash -r "addpath '+matlab_scripts_loc+';imhistmatch_script'
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['result_filename'] = os.path.abspath(self.inputs.result_filename)   
        return outputs 
    
#==============================================================================    

#==============================================================================
# gets the echo times from a procpar file (for Varian data)
#==============================================================================
def get_te_from_procpar_fn(filename):
    import numpy as np
    esp=0
    espincr=0    
    with open(filename) as f:
        #read_data=f.read()
        for line in f:            
            if line.startswith('te '):            
                te=float(next(f).split(" ")[1])            
            elif line.startswith('esp '):            
                esp=float(next(f).split(" ")[1])
            elif line.startswith('espincr '):            
                espincr=float(next(f).split(" ")[1])
            elif line.startswith('ne '):            
                ne=int(next(f).split(" ")[1])
    f.close()
    te_list=te + np.arange(ne)*esp + (1e-6*np.arange(-1,ne-1))* (np.arange(ne)*espincr/2.0)
    te_list=list(te_list)
    return te_list


GetTeFromProcpar = Function(input_names=["filename"],
                             output_names=["te"],
                             function=get_te_from_procpar_fn)
#==============================================================================    

#==============================================================================
# gets center frequency from a procpar file (for Varian data)
#==============================================================================
def get_CF_from_procpar_fn(filename):    
    
    with open(filename) as f:
        #read_data=f.read()
        for line in f:            
            if line.startswith('sfrq '):            
                CF=float(next(f).split(" ")[1])*1e6
    f.close()    
    return CF
GetCFFromProcpar = Function(input_names=["filename"],
                             output_names=["cf"],
                             function=get_CF_from_procpar_fn)
#==============================================================================    

#==============================================================================
#gets center frequency from json files (for BIDS)
#==============================================================================
def get_CF_from_json_fn(filename): 
    
    import json
    if type(filename)==list:
        filename=filename[0]
    with open(filename) as f:
        CF=float(json.load(f)['ImagingFrequency'])*1e6       
    f.close()    
    return CF
GetCFFromJson = Function(input_names=["filename"],
                             output_names=["cf"],
                             function=get_CF_from_json_fn)

#==============================================================================    
 

#==============================================================================
#open list of mag images and return the average and weights based on local
#neighbourhood snr
#==============================================================================
class GetAvgAndWeightsFromMagInputSpec(CommandLineInputSpec):

    mag = InputMultiPath(File(exists=True), desc='A list of magnitude echoes',                                 
                                 copyfile=False, mandatory=True)
    snr_window_sz = traits.Float(15,desc='The size of the window used to calculate local SNR in mm',
                                 mandatory=False,usedefault=True)
    avg_out_filename = File(desc="Output file",mandatory=True)
    weight_out_filename = File(desc="Output file",mandatory=True)
    
class GetAvgAndWeightsFromMagOutputSpec(TraitedSpec):    
    avg_out_filename = File(desc="Output file containing the estimated Frequency",
                   exists=True)
    weight_out_filename = File(desc="Output file containing the estimated Frequency",
                   exists=True)
    
class GetAvgAndWeightsFromMag(BaseInterface):
    input_spec = GetAvgAndWeightsFromMagInputSpec
    output_spec = GetAvgAndWeightsFromMagOutputSpec    
    
    def _run_interface(self, runtime):
        mag_filename_list=self.inputs.mag
        snr_window_sz_mm=self.inputs.snr_window_sz
        avg_out_filename=self._list_outputs()['avg_out_filename']
        weight_out_filename=self._list_outputs()['weight_out_filename']        
        ne=len(mag_filename_list)
        mag_img_obj=nib.load(mag_filename_list[0])
        shape=mag_img_obj.get_shape()
        voxel_size=mag_img_obj.header['pixdim'][1:4]
        snr_window_sz_vxl=np.ceil(snr_window_sz_mm/np.array(voxel_size)).astype('int')        
        avg=np.zeros(shape)
        weight_avg=0
        weight=np.empty(shape+(ne,))
        count=0
        mag_filename_list.sort()
        for imgloc in mag_filename_list:
            curr_img=(nib.load(imgloc).get_data()).astype('float')
            curr_img_avg=np.abs(curr_img).mean()            
            local_mean = ndimage.uniform_filter(curr_img, snr_window_sz_vxl)
            local_sqr_mean = ndimage.uniform_filter(curr_img**2, snr_window_sz_vxl)
            local_var = local_sqr_mean - local_mean**2    
            weight[...,count]=local_mean/local_var #snr
            avg+=curr_img_avg*curr_img #weighted avg, later echos weighted less
            weight_avg+=curr_img_avg
            count+=1
        avg=avg/weight_avg        
        
        niftifile=nib.Nifti1Pair(avg,mag_img_obj.affine)
        nib.save(niftifile,avg_out_filename)
        niftifile=nib.Nifti1Pair(weight,mag_img_obj.affine)
        nib.save(niftifile,weight_out_filename)
        
        return runtime
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['avg_out_filename'] = os.path.abspath(self.inputs.avg_out_filename)
        outputs['weight_out_filename'] = os.path.abspath(self.inputs.weight_out_filename)
        return outputs 
#==============================================================================    

#==============================================================================
#CalculatReliabilityMask Interface uses second difference quality map to
#masked out of fieldmap voxels which are noisey.
#Masking worked better than weighting fieldmap voxel by its reliability
#==============================================================================
class CalculatReliabilityMaskInputSpec(CommandLineInputSpec):
    phase = InputMultiPath(File(exists=True), desc='A list of phase echoes',                                 
                                 copyfile=False, mandatory=True)    
    mask = File(desc = 'Mask (3D) where to calculate phase reliability', 
                mandatory=True)
    threshold = traits.Float(3,desc='Threshold for second difference quality map',
                             mandatory=False,usedefault=True)
    reliability_mask_filename = File(desc="Output file",mandatory=True)   
    reliability_filename = File(desc="Reliability output file",mandatory=True)     
    
class CalculatReliabilityMaskOutputSpec(TraitedSpec):    
    reliability_mask_filename = File(desc="Phase reliability mask",exists=True)
    reliability_filename = File(desc="The reliability map used to make the mask",exists=True)
    
class CalculatReliabilityMask(BaseInterface):
    input_spec = CalculatReliabilityMaskInputSpec
    output_spec = CalculatReliabilityMaskOutputSpec
    
    def _run_interface(self, runtime):
        phase_filename_list=self.inputs.phase        
        mask_filename=self.inputs.mask       
        threshold=self.inputs.threshold
        reliability_mask_filename=self._list_outputs()['reliability_mask_filename']
        reliability_filename=self._list_outputs()['reliability_filename']
        ne=len(phase_filename_list)
        ph_img_obj=nib.load(phase_filename_list[0])        
        shape=ph_img_obj.get_shape()        
        ph_imgs=np.empty(shape+(ne,))        
        for echo in range(ph_imgs.shape[-1]):
            ph_imgs[...,echo]=(nib.load(phase_filename_list[echo]).get_data())            
        reliability=np.empty_like(ph_imgs)
        mask=(nib.load(mask_filename).get_data()).astype('bool')
        for echo in range(ph_imgs.shape[-1]):
            reliability[...,echo]=cr.calculateReliability(ph_imgs[...,echo].astype('float32'),mask)
                
        mask=mask[...,np.newaxis]
        mask.repeat(ne,axis=-1) 
        """
        #absolute value of second difference dependent on image mean
        #dependent on background field
        #could try using std units to avoid this
        threshold_std=self.inputs.threshold_std_units
        tmpmean=reliability[mask].mean()
        tmpmask=(reliability<tmpmean)*(reliability>0)
        mean=reliability[tmpmask].mean()        
        std=reliability[tmpmask].std()
        threshold=mean + threshold_std * std  
        """        
        reliability_mask=(reliability<threshold)*mask
        reliability_mask=np.squeeze(reliability_mask)
        reliability_mask=reliability_mask.astype('int8')

        """        
        #smooth reliability mask to avoid noisey on/off pixel patterns
        rtmp=reliability_mask.copy()
        kernelSz_mm=0.5#mm
        KernelSz_vxl=np.ceil(kernelSz_mm/np.array(voxelSize)).astype('int')        
        for i in range(5):
            rtmp=ndimage.gaussian_filter(rtmp.astype('float'),np.concatenate((KernelSz_vxl,(0,))))
            rtmp=rtmp*reliabilityMask          
        rtmp=ndimage.gaussian_filter(rtmp.astype('float'),(0,0,0,.5))
        #rtmp=rtmp*reliabilityMask
        for echo in range(reliabilityMask.shape[-1]):
            rtmp[...,echo]=ndimage.binary_erosion(rtmp[...,echo])
        reliability_mask=rtmp
        #"""
        
        niftifile=nib.Nifti1Pair(reliability_mask,ph_img_obj.affine)
        nib.save(niftifile,reliability_mask_filename)   
        niftifile=nib.Nifti1Pair(reliability,ph_img_obj.affine)
        nib.save(niftifile,reliability_filename)                       
        return runtime
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['reliability_mask_filename'] = os.path.abspath(self.inputs.reliability_mask_filename)       
        outputs['reliability_filename'] = os.path.abspath(self.inputs.reliability_filename)       
        return outputs 
#==============================================================================

#==============================================================================
"""
The original mask will be eroded in order to constrain where to do reliability 
calculation.  To erode the mask it will be inverted and convolved with spherical 
kernel of user defined size and the thresholded.  Reliability will be calculated
in the contsrained section between original mask and eroded mask - any voxels
that are unreliable will be removed in order to createa  new mask.
"""
#==============================================================================

class TrimMaskUsingReliabilityInputSpec(CommandLineInputSpec):
    phase = File(exists=True, desc='Phase image used to calculate reliability',
                 copyfile=False, mandatory=True)    
    mask = File(desc = 'Mask to be eroded',mandatory=True)
    erosion_sz = traits.Float(5,desc='Erosion in mm',mandatory=False,usedefault=True)
    threshold = traits.Float(0,desc='Threshold for second difference quality map',
                             mandatory=False,usedefault=True)    
    trimmed_mask_filename = File(desc="Mask output file",mandatory=True)    
    reliability_filename = File(desc="Reliability output file",mandatory=True) 
    
class TrimMaskUsingReliabilityOutputSpec(TraitedSpec):    
    trimmed_mask_filename = File(desc="The trimmed mask",exists=True)
    reliability_filename = File(desc="The reliability map used to make the mask",exists=True)
    
class TrimMaskUsingReliability(BaseInterface):
    input_spec = TrimMaskUsingReliabilityInputSpec
    output_spec = TrimMaskUsingReliabilityOutputSpec
    
    def _run_interface(self, runtime):
        phaseFilename=self.inputs.phase        
        maskFilename=self.inputs.mask       
        threshold=self.inputs.threshold
        rad=self.inputs.erosion_sz
        trimmed_mask_filename=self._list_outputs()['trimmed_mask_filename']
        reliability_filename=self._list_outputs()['reliability_filename']
                
        ph_img_obj=nib.load(phaseFilename)
        voxel_size=ph_img_obj.header['pixdim'][1:4]        
        ph_img=ph_img_obj.get_data()               
        mask=nib.load(maskFilename).get_data()
        
        px,py,pz=int(rad/voxel_size[0]),int(rad/voxel_size[1]),int(rad/voxel_size[2])
        x=np.linspace(-rad,rad,2*px+1)
        y=np.linspace(-rad,rad,2*py+1)
        z=np.linspace(-rad,rad,2*pz+1)
        Y,X,Z=np.meshgrid(y,x,z)
        dist2=(X**2+Y**2+Z**2)
        circ=dist2<=rad**2        
        
        mask_eroded=ndimage.binary_erosion(mask,circ)      
        mask_dilated=ndimage.binary_dilation(mask,iterations=3)
        reliability=cr.calculateReliability(ph_img.astype('float32'),mask_dilated.astype('bool'))
        
        """
        #absolute value of second difference dependent on image mean
        #dependent on background field
        #could try using std units to avoid this
        threshold_std=self.inputs.threshold_std
        mean=reliability[mask.astype('bool')].mean()
        std=reliability[mask.astype('bool')].std()        
        threshold=mean + threshold_std * std
        threshold=threshold_std
        """        
                
        newmask=(((reliability<threshold)*mask+mask_eroded)>0)
        newmask=ndimage.binary_fill_holes(newmask)
        newmask=ndimage.binary_opening(newmask).astype('int8')        
      
        niftifile=nib.Nifti1Pair(newmask,ph_img_obj.affine)
        nib.save(niftifile,trimmed_mask_filename)       
        niftifile=nib.Nifti1Pair(reliability,ph_img_obj.affine)
        nib.save(niftifile,reliability_filename)                
        return runtime
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['trimmed_mask_filename'] = os.path.abspath(self.inputs.trimmed_mask_filename)       
        outputs['reliability_filename'] = os.path.abspath(self.inputs.reliability_filename)       
        return outputs 
#==============================================================================
    
#==============================================================================
#use pyQSM package to estimate the frequency from phase images
#==============================================================================
class EstimateFrequncyFromWrappedPhaseInputSpec(CommandLineInputSpec):

    phase = InputMultiPath(File(exists=True), desc='A list of phase echoes',
                                 argstr="--phase %s", position=1,
                                 copyfile=False, mandatory=True)
    json = InputMultiPath(File(exists=True), desc='A list of json files containing the tes',
                                 argstr="--phase %s", position=1,
                                 copyfile=False, mandatory=True)    
    truncate_echo = traits.Int(0,desc='Cut off echoes after this point',argstr='-truncate %s',
                          position=3,mandatory=False,usedefault=True)    
    mask = File(desc = 'Mask (3D) where to estimate frequency', argstr="--mask %s",
                  position=4,mandatory=True)    
    weight = File(desc = 'Wegiths (4D) for phase measurements', argstr="--weight %s",
                  position=5,mandatory=False)
    
    freq_filename = File(desc="Output file",argstr='--freq %s',position=6,               
                   mandatory=True)
    
class EstimateFrequncyFromWrappedPhaseOutputSpec(TraitedSpec):    
    freq_filename = File(desc="Output file containing the estimated Frequency",
                   exists=True)
    
class EstimateFrequncyFromWrappedPhase(BaseInterface):
    input_spec = EstimateFrequncyFromWrappedPhaseInputSpec
    output_spec = EstimateFrequncyFromWrappedPhaseOutputSpec
    
    def _run_interface(self, runtime):        
        phase_filename_list=self.inputs.phase
        json_filename_list=self.inputs.json
        truncate_echo=self.inputs.truncate_echo
        if not truncate_echo:
            truncate_echo=None
        mask_filename=self.inputs.mask       
        freq_filename=self._list_outputs()['freq_filename']       
               
        if isdefined(self.inputs.weight):
            weight_filename=self.inputs.weight
            weight=nib.load(weight_filename).get_data()
        else:
            weight=None
        ne=len(phase_filename_list)
        ph_img_obj=nib.load(phase_filename_list[0])
        shape=ph_img_obj.get_shape()
        voxel_size=ph_img_obj.header['pixdim'][1:4]
        ph_img=np.empty(shape+(ne,))
        if weight is None:
            weight=np.ones(shape+(ne,))
        count=0
        phase_filename_list.sort()
        json_filename_list.sort()
        te=np.empty(ne)
        for imgloc,jsonloc in zip(phase_filename_list,json_filename_list):            
            #assert os.path.splitext(jsonloc)[0] in os.path.splitext(imgloc)[0]
            ph_img[...,count]=nib.load(imgloc).get_data()
            with open(jsonloc) as f:
                te[count]=float(json.load(f)['EchoTime'])          
            count+=1               
        mask=nib.load(mask_filename).get_data()
        freq=pyQSM.frequencyEstimate.estimateFrequencyFromWrappedPhase(ph_img,voxel_size,te,mask,weight,truncateEcho=truncate_echo)                
        niftifile=nib.Nifti1Pair(freq,ph_img_obj.affine)
        nib.save(niftifile,freq_filename)           
        return runtime
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['freq_filename'] = os.path.abspath(self.inputs.freq_filename)        
        return outputs 
#==============================================================================    

#==============================================================================
# Cornell implementation of resharp, executed in matlab
#==============================================================================
class RESHARPInputSpec(CommandLineInputSpec):

    freq = File(exists=True, desc='Input filename',
                                 argstr='%s', position=1,
                                 copyfile=False, mandatory=True)        
    mask = File(exists=True, desc = 'Mask (3D)', argstr='%s',
                  position=2,copyfile=False,mandatory=True)
    
    radius = traits.Float(5.0,desc='the radius of the spherical mean value operation in mm',argstr='%s',
                          position=4,mandatory=False,usedefault=True)
    
    alpha = traits.Float(0.05,desc='the regularizaiton parameter used in Tikhonov',argstr='%s',
                          position=5,mandatory=False,usedefault=True)    
    
    LFS_filename = File(desc="Output LFS filename",argstr='%s',position=6,               
                   mandatory=True)   
    
    quit_matlab = traits.String(';quit;"',desc='needed to quit matlab',argstr='%s',
                      position=7,mandatory=False,usedefault=True)   
    
class RESHARPOutputSpec(TraitedSpec):    
    LFS_filename = File(desc="Local frequency shift with background removed",
                   exists=True)
    LFS_mask_filename = File(desc="Local frequency shift with background removed",
                   exists=True)
    
class RESHARP(CommandLine):
    input_spec = RESHARPInputSpec
    output_spec = RESHARPOutputSpec
    cmd = 'matlab -nodisplay -nosplash -r "addpath '+matlab_scripts_loc+';RESHARP_script'
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['LFS_filename'] = os.path.abspath(self.inputs.LFS_filename)
        LFS_mask_filename=outputs['LFS_filename']
        LFS_mask_filename=os.path.splitext(LFS_mask_filename)[0]
        LFS_mask_filename=LFS_mask_filename+"_mask.nii"
        outputs['LFS_mask_filename'] = os.path.abspath(LFS_mask_filename)
        return outputs 
#==============================================================================

#==============================================================================
# Western University CFMM implementation of MEDI, pyQSM package
#==============================================================================

class MEDI_CFMM_implementation_InputSpec(CommandLineInputSpec):

    lfs_loc = File(exists=True, desc='Input lfs/rdf filename',
                                 argstr="%s", position=1,
                                 copyfile=False, mandatory=True)
    iMag_loc = File(exists=True, desc='Input avg magnitude for edge finding',
                                 argstr="%s", position=2,
                                 copyfile=False, mandatory=True)
    weight_loc = File(exists=True, desc='Weights for data fidelity in MEDI',
                                 argstr="%s", position=6,
                                 copyfile=False, mandatory=False)
    mask_loc = File(exists=True, desc = 'Mask (3D)', argstr="%s",
                  position=3,copyfile=False,mandatory=True)
    lamda = traits.Float(2.0,desc='Regularization parameter for MEDI',argstr='%s',
                          position=4,mandatory=False,usedefault=True)
    #note that MEDI wants it's rdf as lfs*te where te~0.004.  
    #We use data that is not multiplied by te and so our data is 3 orders of magnitude larger
    #therefore if in matlab MEDI you use lambda of 500, then in this implementation try using 500*.004~2
    maxiter =  traits.Int(200,desc='Max number of interations to perform in optimization',argstr='%s',
                          position=3,mandatory=False,usedefault=True)
    tol =  traits.Float(0.0015,desc='Tolerance. Quit optimization loop when change in susceptibility estimate is less than tolerance.',argstr='%s',
                          position=3,mandatory=False,usedefault=True)
    CF = traits.Float(298060000.0,desc='Center frequency, used to return result in ppb',argstr='%s',
                          position=4,mandatory=False,usedefault=True)
    susceptibility_filename = File(desc="Output susceptibility filename",argstr='%s',position=5,               
                   mandatory=True)    
   
    
class MEDI_CFMM_implementation_OutputSpec(TraitedSpec):    
    susceptibility_filename = File(desc="Susceptibility estimate",exists=True)
    
class MEDI_CFMM_implementation(BaseInterface):
    input_spec = MEDI_CFMM_implementation_InputSpec
    output_spec = MEDI_CFMM_implementation_OutputSpec
    def _run_interface(self, runtime):
        
        rdf_file=self.inputs.lfs_loc
        mask_file=self.inputs.mask_loc
        iMag_file=self.inputs.iMag_loc
        WData_file=self.inputs.weight_loc
        alpha=self.inputs.lamda #split bregman alpha is the weight on the data fidelity, MEDI uses lambda for this
        lamda=2*alpha #split bregman uses lambda to refer to the weight which keeps bregman parameters close to their corresponding data
        maxiter=self.inputs.maxiter
        tol=self.inputs.tol
        CF=self.inputs.CF
        
        rdf_obj=nib.load(rdf_file)
        rdf=rdf_obj.get_data()
        voxel_size=rdf_obj.header['pixdim'][1:4]
        mask=nib.load(mask_file).get_data()
        iMag=nib.load(iMag_file).get_data().transpose(1,0,2)[:,::-1,:]
        WData_sqrd=nib.load(WData_file).get_data()[...,0]**2*mask
        
        #create mask to protect edges from regularization
        percentage=0.1 #keep 10% of voxels with highest gradient value
        WGradx,WGrady,WGradz=np.gradient(iMag,*voxel_size)
        wG=np.sqrt(WGradx**2+WGrady**2+WGradz**2)
        thresh=wG[mask.astype('bool')].max()
        num=((mask*wG)>thresh).sum()
        den=(mask>0).sum()        
        while float(num)/den<percentage:    
            thresh=thresh*0.95
            num=((mask*wG)>thresh).sum()
            den=(mask>0).sum()
        wG=(mask*wG)<thresh
        WGradx=WGrady=WGradz=wG*mask
        
        #use CFMM local implementation of MEDI
        #when using matlab implementation of MEDI solver, weighting matrix causes numerical artifacts
        FOV=list(voxel_size*rdf.shape)
        Xmap=pyQSM.Liu2012.SplitBregmanL1Grad(rdf,FOV, WData_sqrd, WGradx,WGrady,WGradz,alpha,lamda,maxiter=maxiter,tol=tol)
        Xmap=Xmap/((CF*1e-9)*2*np.pi)*mask #part per billion

        susceptibility_filename=self._list_outputs()['susceptibility_filename']               
        niftifile=nib.Nifti1Pair(Xmap,rdf_obj.affine)
        nib.save(niftifile,susceptibility_filename)        
        return runtime
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['susceptibility_filename'] = os.path.abspath(self.inputs.susceptibility_filename)        
        return outputs 
#==============================================================================

#==============================================================================
#Cornell implementation of MEDI, implemented in matlab
#==============================================================================

class MEDIInputSpec(CommandLineInputSpec):

    lfs_loc = File(exists=True, desc='Input LFS filename',
                                 argstr="%s", position=1,
                                 copyfile=False, mandatory=True)#, xor=['complex_names'])  
    iMag_loc = File(exists=True, desc='Input avg magnitude for edge finding',
                                 argstr="%s", position=2,
                                 copyfile=False, mandatory=True)#, xor=['complex_names'])  
    weight_loc = File(exists=True, desc='Weights for data fidelity in MEDI',
                                 argstr="%s", position=6,
                                 copyfile=False, mandatory=False)#, xor=['complex_names'])        
    mask_loc = File(exists=True, desc = 'Mask (3D)', argstr="%s",
                  position=3,copyfile=False,mandatory=True)
    lamda = traits.Float(200.0,desc='Regularization parameter for MEDI',argstr='%s',
                          position=4,mandatory=False,usedefault=True)    
    susceptibility_filename = File(desc="Output susceptibility filename",argstr='%s',position=5,               
                   mandatory=True)    
    quit_matlab = traits.String(';quit;"',desc='needed to quit matlab',argstr='%s',
                      position=7,mandatory=False,usedefault=True)   
    
class MEDIOutputSpec(TraitedSpec):    
    susceptibility_filename = File(desc="Susceptibility estimate",exists=True)
    
class MEDI(CommandLine):
    input_spec = MEDIInputSpec
    output_spec = MEDIOutputSpec
    cmd = 'matlab -nodisplay -nosplash -r "addpath '+matlab_scripts_loc+';MEDI_script'
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['susceptibility_filename'] = os.path.abspath(self.inputs.susceptibility_filename)        
        return outputs 
#==============================================================================
    
#==============================================================================
#runs a matlab script that runs SS_TV for a range of lagrange parameter values.
#uses elbow method to find optimal lagrange parameter. modified from bilgic code
#==============================================================================
class SBFindOptParamInputSpec(CommandLineInputSpec):
    freq_loc = File(exists=True, desc='Input freq filename',
                                 argstr="%s", position=1,
                                 copyfile=False, mandatory=True)    
    reliability_mask_loc = File(exists=True, desc='Phase voxels to leave out of data fidelity term',
                                 argstr="%s", position=3,
                                 copyfile=False, mandatory=False) 
    mask_loc = File(exists=True, desc = 'Mask (3D)', argstr="%s",
                  position=2,copyfile=False,mandatory=True)
    result_mat_file_loc = File(desc="Output mat filename where optimal param is stored",argstr='%s',position=4,               
                   mandatory=True)    
    quit_matlab = traits.String(';quit;"',desc='needed to quit matlab',argstr='%s',
                      position=7,mandatory=False,usedefault=True)   
    
class SBFindOptParamOutputSpec(TraitedSpec):    
    result_mat_file_loc = File(desc=".mat file storing optimal SB parameter",exists=True)
    
class SBFindOptParam(CommandLine):
    input_spec = SBFindOptParamInputSpec
    output_spec = SBFindOptParamOutputSpec
    cmd = 'matlab -nodisplay -nosplash -r "addpath '+matlab_scripts_loc+';LCurveOptimalParam_script'
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['result_mat_file_loc'] = os.path.abspath(self.inputs.result_mat_file_loc)        
        return outputs 
#==============================================================================

    
#==============================================================================
#runs Bilgic lab's matlab implementation of SS_TV 
#==============================================================================
class SS_TVInputSpec(CommandLineInputSpec):
    freq_loc = File(exists=True, desc='Input freq filename',
                                 argstr="%s", position=1,
                                 copyfile=False, mandatory=True)    
    reliability_mask_loc = File(exists=True, desc='Phase voxels to leave out of data fidelity term',
                                 argstr="%s", position=3,
                                 copyfile=False, mandatory=False) 
    mask_loc = File(exists=True, desc = 'Mask (3D)', argstr="%s",
                  position=2,copyfile=False,mandatory=True)
    CF = traits.Float(298060000.0,desc='Center frequency, used to return result in ppb',argstr='%s',
                          position=4,mandatory=False,usedefault=True)
    alpha = traits.Float(0.3,desc='Regularization parameter for MEDI',argstr='%s',
                          position=5,mandatory=False,usedefault=True)
    B0_dir = traits.Int(3,desc='B0 direction (1, 2, or 3)',argstr='%s',
                          position=6,mandatory=False,usedefault=True)
    susceptibility_filename = File(desc="Output susceptibility filename",argstr='%s',position=7,               
                   mandatory=True)    
    quit_matlab = traits.String(';quit;"',desc='needed to quit matlab',argstr='%s',
                      position=8,mandatory=False,usedefault=True)   
    
class SS_TVOutputSpec(TraitedSpec):    
    susceptibility_filename = File(desc="Susceptibility estimate",exists=True)
    
class SS_TV(CommandLine):
    input_spec = SS_TVInputSpec
    output_spec = SS_TVOutputSpec
    cmd = 'matlab -nodisplay -nosplash -r "addpath '+matlab_scripts_loc+';SS_TV_script'
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['susceptibility_filename'] = os.path.abspath(self.inputs.susceptibility_filename)        
        return outputs     
    
class SS_TV_mcr(CommandLine):
    input_spec = SS_TVInputSpec
    output_spec = SS_TVOutputSpec
    _cmd = matlab_scripts_loc+'/run_SS_TV_script.sh '+mcr_location+' '
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['susceptibility_filename'] = os.path.abspath(self.inputs.susceptibility_filename)        
        return outputs        
#==============================================================================

#==============================================================================
#useful for creating filenames from paths
#==============================================================================
def replace_slash_fn(filename):  
        renamed="_".join(str(filename).split("/"))            
        return renamed
replace_slash = Function(input_names=["filename"],
                             output_names=["renamed"],
                             function=replace_slash_fn)
#==============================================================================
  
    
if __name__=='__main__':
    pass
    

    
