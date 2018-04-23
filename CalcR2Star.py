#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 13:39:05 2017

@author: akuurstr
"""
import multiprocessing as mp
from multiprocessing import Pool,sharedctypes
import numpy as np
from scipy.optimize import leastsq
import nibabel as nib

def model_monoexponential_R2Star(TE,A,B,R2star,f):
    #do not use mag/phase for scaling parameter A!!! use re/imag
    return (A+1j*B)*np.exp(-TE*R2star+1j*f*TE)
    
def residuals(params,TE,data):
    diff=model_monoexponential_R2Star(TE,*params)-data
    diffconcat=np.empty(diff.size*2,dtype='float')
    diffconcat[0:diffconcat.size:2]= diff.real
    diffconcat[1:diffconcat.size:2]= diff.imag   
    return diffconcat

def calcR2star_parallelizable_unit(indx):
    x0=(initial_r2_guess_scale_real[indx],initial_r2_guess_scale_imag[indx],initial_r2_guess_exponent_real[indx],initial_r2_guess_exponent_imag[indx])    
    complex_signal=[]
    for i in range(len(TE)):                
        exec("complex_signal.append(echo%s_real[%s]+1j*echo%s_imag[%s])" % (i,indx,i,indx))
 
    complex_signal=np.array(complex_signal)
    
    lstsqopt=leastsq(residuals,x0,args=(np.array(TE[:]),complex_signal),maxfev=10000,full_output=True)
    
    ss_err=(lstsqopt[2]['fvec']**2).sum()
    ss_tot=(np.abs(complex_signal-complex_signal.mean())**2).sum()
    rsquared=1-(ss_err/ss_tot)  
    
    if lstsqopt[-1] in (1,2,3,4):
        return lstsqopt[0][2],rsquared
    else:        
        return np.nan,np.inf

    
def _init_R2star_shared_arrays(real_echoes, imag_echoes, initial_r2_guess_scale_real_shared,initial_r2_guess_scale_imag_shared, initial_r2_guess_exponent_real_shared,initial_r2_guess_exponent_imag_shared, TE_shared):
    """ Each pool process calls this initializer. Load the array to be populated into that process's global namespace """
    
    for i in range(len(TE_shared)):
            exec("global echo%s_real; echo%s_real=real_echoes[i]" % (i,i))
            exec("global echo%s_imag; echo%s_imag=imag_echoes[i]" % (i,i))
                
    global initial_r2_guess_scale_real
    initial_r2_guess_scale_real=initial_r2_guess_scale_real_shared    
    global initial_r2_guess_scale_imag
    initial_r2_guess_scale_imag=initial_r2_guess_scale_imag_shared
    global initial_r2_guess_exponent_real
    initial_r2_guess_exponent_real=initial_r2_guess_exponent_real_shared    
    global initial_r2_guess_exponent_imag
    initial_r2_guess_exponent_imag=initial_r2_guess_exponent_imag_shared
    
    global TE
    TE = TE_shared



    
def calc_R2star_fn(img,mask,TE):
    img_shape=img.shape
    img=img/(img[...,0]*mask).max() #useful if you want to use ftol to exit the curve_fit
    img=img.reshape(np.prod(img_shape[:-1]),img_shape[-1])    

    A=np.array([np.ones_like(TE),-TE]).T   
    y=np.log(img.reshape(np.prod(img.shape[:-1]),img.shape[-1])+1e-20).T 
    x=np.dot(np.linalg.pinv(A),y)    
    
    indices=np.arange(np.prod(img_shape[:-1]))
    indices=indices.reshape(img_shape[:-1])    
    indices=indices[mask.astype('bool')]
        
    real_echoes=[]
    imag_echoes=[]
    
    for i in range(len(TE)):
            exec("echo%s_real_shared=sharedctypes.Array('d',img[...,i].real.ravel(), lock=False);real_echoes.append(echo%s_real_shared)" % (i,i))
            exec("echo%s_imag_shared=sharedctypes.Array('d',img[...,i].imag.ravel(), lock=False);imag_echoes.append(echo%s_imag_shared)" % (i,i))
    
    x_scale_real_shared=sharedctypes.Array('d',x[0,...].real.ravel(),lock=False)        
    x_scale_imag_shared=sharedctypes.Array('d',x[0,...].imag.ravel(),lock=False)    
    x_exponent_real_shared=sharedctypes.Array('d',x[1,...].real.ravel(),lock=False)        
    x_exponent_imag_shared=sharedctypes.Array('d',x[1,...].imag.ravel(),lock=False)    
    
    TE_shared=sharedctypes.Array('d',TE,lock=False)    
    
    mp.set_start_method('spawn')   
    chunksize=10000
    maxtasksperchild=5
    
    pool = Pool(maxtasksperchild=maxtasksperchild, initializer=_init_R2star_shared_arrays,initargs=(real_echoes, imag_echoes,x_scale_real_shared,x_scale_imag_shared, x_exponent_real_shared,x_exponent_imag_shared,TE_shared))        
    result=np.array(pool.map(calcR2star_parallelizable_unit,indices,chunksize=chunksize))
        
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


if __name__ == '__main__':
    from argparse import ArgumentParser 
    import json    
    import os
    parser = ArgumentParser()
    parser.add_argument('--mag_filename_list',required=True,nargs="+")
    parser.add_argument('--phase_filename_list',required=True,nargs="+")
    parser.add_argument('--json_filename_list',required=True,nargs="+")
    parser.add_argument('--freq_filename',required=True)
    parser.add_argument('--mask_filename',required=True)
    parser.add_argument('--r2star_filename',required=True)
    parser.add_argument('--neg_mask_filename',required=True)
    parser.add_argument('--nan_mask_filename',required=True)
    args = parser.parse_args()
    
    #inputs from command line
    mag_filename_list=args.mag_filename_list
    phase_filename_list=args.phase_filename_list
    json_filename_list=args.json_filename_list
    freq_filename=args.freq_filename
    mask_filename=args.mask_filename
         
    #output names from command line
    r2star_filename=args.r2star_filename
    neg_mask_filename=args.neg_mask_filename
    nan_mask_filename=args.nan_mask_filename
               
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
    data_reduced_phase=complex_img*np.exp(-1j*freq[...,np.newaxis]*te)     
    R2star_img,goodness_of_fit,neg_mask,nan_mask=calc_R2star_fn(data_reduced_phase,mask,te)    
    neg_mask=neg_mask.astype('int8')
    nan_mask=nan_mask.astype('int8')
    
    niftifile=nib.Nifti1Pair(R2star_img,mag_img_obj.affine)
    nib.save(niftifile,r2star_filename)
    niftifile=nib.Nifti1Pair(goodness_of_fit,mag_img_obj.affine)
    #extension should start at leading dot instead of last dot
    #use splitext in reverse
    ext,name=os.path.splitext(r2star_filename[::-1])
    name=name[:0:-1]
    ext='.'+ext[::-1]    
    nib.save(niftifile,name+'_fit'+ext)
    niftifile=nib.Nifti1Pair(neg_mask,mag_img_obj.affine)
    nib.save(niftifile,neg_mask_filename)
    niftifile=nib.Nifti1Pair(nan_mask,mag_img_obj.affine)
    nib.save(niftifile,nan_mask_filename)
        
 
  

    
    
    