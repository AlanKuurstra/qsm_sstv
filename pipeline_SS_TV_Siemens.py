#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 11:07:39 2017

@author: akuurstr

-in custom python nodes, using nibabel to open niftis and always copy the input
 affine matrix to the output nifti.
-in matlab nodes, using load_untouch_nii to load data without applying
 affine transformations. output nifti is saved using input's header
-user is expected to provide which dimension the nifti file stores the 
 B0 direction
"""
import numpy as np
from nipype.interfaces import utility as niu
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.fsl as fsl
from interfaces_CFMM import GetAvgAndWeightsFromMag,CalculatReliabilityMask,\
TrimMaskUsingReliability,EstimateFrequncyFromWrappedPhase,GetCFFromJson,SS_TV,SS_TV_mcr,\
CalcR2Star,CalcR2Star_cmd,SiemensPhasePreprocess,replace_slash
from nipype.interfaces.utility import Function
from bids.grabbids import BIDSLayout
import os

def create_pipeline_SS_TV(  bids_dir,
                            work_dir,
                            out_dir,
                            subjects,
                            sessions,
                            mag_match_pattern,
                            phase_match_pattern,
                            keep_unnecessary_outputs,
                            FAST_bias_iters,
                            FAST_bias_lowpass,
                            FAST_num_classes,
                            BET_frac,
                            freq_weights__snr_window_sz,
                            truncate_echo,
                            SS_TV_lagrange_parameter,
                            B0_dir,
                            scnd_diff_reliability_thresh_trim,
                            scnd_diff_reliability_thresh_noise):
    layout=BIDSLayout(bids_dir)
        
    #can we do this more elegantly?
    first_echo_files=[]    
    for subject in subjects:
        if layout.get_sessions(subject=subject)==[]:
            if sessions==['.*']:
                first_echo_files=first_echo_files+layout.get(subject=subject,modality='anat',extensions='.*part-phase.*echo-0*1.*.nii.*',)
            else:
                print("Warning: Session filter applied, but subject "+subject+" has no bids session information. This subject has been ignored.")
        else:
            for session in sessions:
                first_echo_files=first_echo_files+layout.get(subject=subject,session=session,modality='anat',extensions='.*part-phase.*echo-0*1.*.nii.*',)        
    anat_folders=[]
    for img in first_echo_files:
        full_dirname=os.path.dirname(img.filename)
        remove_base_dir=full_dirname.replace(bids_dir,'')
        remove_leading_slash=remove_base_dir.lstrip(os.sep)
        anat_folders.append(remove_leading_slash)
    list(set(anat_folders)).sort()
        
    #IdentityInterface is useful for passing subject directory structure to datasink
    infosource = pe.Node(niu.IdentityInterface(fields=['subject_id']), name="infosource")
    infosource.iterables = ('subject_id', anat_folders)
    
    ### NODES AND PARAMETERS
    datasource=pe.Node(nio.DataGrabber(infields=['subject_id'],outfields=['phase_images','mag_images','phase_jsons','mag_jsons']),name='datasource')
    datasource.inputs.field_template=dict(
            phase_images='%s/'+phase_match_pattern+'.nii*',
            phase_jsons='%s/'+phase_match_pattern+'.json',
            mag_images='%s/'+mag_match_pattern+'.nii*',
            mag_jsons='%s/'+mag_match_pattern+'.json',
            )    
    datasource.inputs.sort_filelist=True
    datasource.inputs.template="*"
    datasource.inputs.base_directory=bids_dir
    
    #this node must change depending on the scanner vendor
    susc_phase_preprocess=pe.Node(SiemensPhasePreprocess(),name='susc_phase_preprocess')
    
    avg_and_freq_estimate_weights=pe.Node(GetAvgAndWeightsFromMag(),name='avg_and_freq_estimate_weights')
    avg_and_freq_estimate_weights.inputs.snr_window_sz=freq_weights__snr_window_sz
    avg_and_freq_estimate_weights.inputs.avg_out_filename="avg.nii.gz"
    avg_and_freq_estimate_weights.inputs.weight_out_filename="weights.nii.gz"          
    
    """
    #spm worked better for varian 7T data
    #if using spm, these prameters are needed
    bias_regularization=.001
    sampling_distance=2.0
    bias_fwhm=30
    
    nonuniformityCorrect_spm=pe.Node(spm.preprocess.Segment(),name='nonuniformityCorrect_spm')
    nonuniformityCorrect_spm.inputs.bias_regularization=bias_regularization
    nonuniformityCorrect_spm.inputs.sampling_distance=sampling_distance
    nonuniformityCorrect_spm.inputs.bias_fwhm=bias_fwhm
    nonuniformityCorrect_spm.inputs.save_bias_corrected=True
    """
    
    nonuniformity_correct_fsl=pe.Node(fsl.FAST(),name='nonuniformity_correct_fsl')       
    nonuniformity_correct_fsl.inputs.img_type=2 #1 for t1, 2 for t2
    nonuniformity_correct_fsl.inputs.bias_iters=FAST_bias_iters #higher for larger nonuniformity
    nonuniformity_correct_fsl.inputs.bias_lowpass=FAST_bias_lowpass #spm uses 30
    nonuniformity_correct_fsl.inputs.number_classes=FAST_num_classes #spm uses 5
    nonuniformity_correct_fsl.inputs.output_biasfield=True
    nonuniformity_correct_fsl.inputs.output_biascorrected=True
    nonuniformity_correct_fsl.interface.estimated_memory_gb = 10
    
    brain_extract=pe.Node(fsl.BET(),name='brain_extract')
    brain_extract.inputs.frac=BET_frac
    brain_extract.inputs.mask=True
    brain_extract.inputs.robust=True
    
    freq_est=pe.Node(EstimateFrequncyFromWrappedPhase(),'freq_est')
    freq_est.inputs.truncate_echo = truncate_echo
    freq_est.inputs.freq_filename = "freq_est.nii.gz"
    freq_est.interface.estimated_memory_gb = 4
    
    R2Star=pe.Node(CalcR2Star_cmd(),'R2Star')
    R2Star.inputs.R2star='R2star.nii.gz'
    R2Star.inputs.neg_mask='negMask.nii.gz'
    R2Star.inputs.nan_mask='nanMask.nii.gz'
    #R2Star.interface.estimated_memory_gb = 5
    
    trim_mask=pe.Node(TrimMaskUsingReliability(),name='trim_mask')
    trim_mask.inputs.erosion_sz = 15.0 #in mm
    trim_mask.inputs.threshold = scnd_diff_reliability_thresh_trim
    trim_mask.inputs.trimmed_mask_filename = "trim_mask.nii.gz"
    trim_mask.inputs.reliability_filename="unreliableMap.nii.gz"
    trim_mask.interface.estimated_memory_gb = 25
    
    unreliable_fieldmap_voxels=pe.Node(CalculatReliabilityMask(),name='unreliable_fieldmap_voxels')
    unreliable_fieldmap_voxels.inputs.threshold=scnd_diff_reliability_thresh_noise
    unreliable_fieldmap_voxels.inputs.reliability_mask_filename="unreliableMask.nii.gz" 
    unreliable_fieldmap_voxels.inputs.reliability_filename="unreliableMap.nii.gz"
    
    CF_value=pe.Node(GetCFFromJson,name='CFValue')   
            
    susceptibility=pe.Node(SS_TV_mcr(),name='susceptibility')
    susceptibility.inputs.quit_matlab='' #use this line when using mcr, comment when using matlab
    susceptibility.inputs.alpha=SS_TV_lagrange_parameter
    susceptibility.inputs.B0_dir=B0_dir
    susceptibility.inputs.susceptibility_filename='susceptibilityMap.nii.gz'
    susceptibility.interface.estimated_memory_gb = 10
    
    fieldmap_reorient=pe.Node(fsl.Reorient2Std(),name='fieldmap_reorient')    
    QSM_reorient=pe.Node(fsl.Reorient2Std(),name='QSM_reorient')
    QSM_brain_mask_reorient=pe.Node(fsl.Reorient2Std(),name='QSM_brain_mask_reorient')
    QSM_noise_mask_reorient=pe.Node(fsl.Reorient2Std(),name='QSM_noise_mask_reorient')
    R2star_reorient=pe.Node(fsl.Reorient2Std(),name='R2star_reorient')
    R2star_fit_reorient=pe.Node(fsl.Reorient2Std(),name='R2star_fit_reorient')
    R2star_neg_mask_reorient=pe.Node(fsl.Reorient2Std(),name='R2star_neg_mask_reorient')   
    
    datasink = pe.Node(nio.DataSink(), name="datasink")    
    datasink.inputs.base_directory = out_dir +'/qsm_sstv/'
    datasink.inputs.parameterization=False    
    
    rename_infosource=pe.Node(replace_slash, "rename_infosource")
    rename_fieldmap = pe.Node(niu.Rename(format_string="%(subject_id)s-fieldmap", keep_ext = True), "rename_fieldmap")
    rename_QSM = pe.Node(niu.Rename(format_string="%(subject_id)s-QSM", keep_ext = True), "rename_QSM")   
    rename_QSM_brain_mask = pe.Node(niu.Rename(format_string="%(subject_id)s-QSM_brainMask", keep_ext = True), "rename_QSM_brain_mask")
    rename_QSM_noise_mask = pe.Node(niu.Rename(format_string="%(subject_id)s-QSM_noiseMask", keep_ext = True), "rename_QSM_noise_mask")
    
    rename_R2star = pe.Node(niu.Rename(format_string="%(subject_id)s-R2star", keep_ext = True), "rename_R2star")   
    rename_R2star_fit = pe.Node(niu.Rename(format_string="%(subject_id)s-R2star_fit", keep_ext = True), "rename_R2star_fit")   
    rename_R2star_neg_mask = pe.Node(niu.Rename(format_string="%(subject_id)s-R2star_negMask", keep_ext = True), "rename_R2star_neg_mask")    
    
    ### PIPELINE CONNECTION
    pipelineDir=work_dir
    wf=pe.Workflow(name="SS_TV")
    wf.base_dir=pipelineDir
    wf.config['execution']['remove_unnecessary_outputs']=False #useful for debugging
    wf.connect([
            (infosource, datasource, [('subject_id', 'subject_id')]),            
            (datasource, avg_and_freq_estimate_weights, [('mag_images', 'mag')]),            
            (datasource, susc_phase_preprocess, [('phase_images', 'infiles')]),
            #spm requires matlab
            #(avg_and_freq_estimate_weights, nonuniformityCorrect_spm, [('avgOutFilename', 'data')]),        
            #(nonuniformityCorrect_spm, brain_extract, [('bias_corrected_image', 'in_file')]),
            (avg_and_freq_estimate_weights, nonuniformity_correct_fsl, [('avg_out_filename', 'in_files')]),        
            (nonuniformity_correct_fsl, brain_extract, [('restored_image', 'in_file')]),
            (susc_phase_preprocess, freq_est, [('outfiles', 'phase')]),
            (datasource, freq_est, [('phase_jsons', 'json')]),
            (brain_extract, freq_est, [('mask_file', 'mask')]),
            (avg_and_freq_estimate_weights, freq_est, [('weight_out_filename', 'weight')]),            
            (freq_est, trim_mask, [('freq_filename', 'phase')]),                                    
            (datasource, R2Star, [('mag_images', 'mag')]),
            (susc_phase_preprocess, R2Star, [('outfiles', 'phase')]),
            (freq_est, R2Star, [('freq_filename','freq_loc')]),            
            (trim_mask, R2Star, [('trimmed_mask_filename','mask')]),
            (datasource, R2Star, [('mag_jsons', 'json')]),        
            (brain_extract, trim_mask, [('mask_file', 'mask')]),
            (freq_est, unreliable_fieldmap_voxels, [('freq_filename','phase')]),            
            (brain_extract, unreliable_fieldmap_voxels, [('mask_file','mask')]),
            (freq_est, susceptibility, [('freq_filename','freq_loc')]), 
            (datasource, CF_value, [('mag_jsons','filename')]),
            (unreliable_fieldmap_voxels, susceptibility, [('reliability_mask_filename','reliability_mask_loc')]),        
            (trim_mask, susceptibility, [('trimmed_mask_filename','mask_loc')]),            
            (CF_value, susceptibility, [('cf','CF')]),
            
            (freq_est, fieldmap_reorient, [('freq_filename', 'in_file')]),
            (susceptibility, QSM_reorient, [('susceptibility_filename', 'in_file')]),
            (trim_mask, QSM_brain_mask_reorient, [('trimmed_mask_filename', 'in_file')]),
            (unreliable_fieldmap_voxels, QSM_noise_mask_reorient, [('reliability_mask_filename', 'in_file')]),
            (R2Star, R2star_reorient, [('R2star', 'in_file')]),
            (R2Star, R2star_fit_reorient, [('R2star_fit', 'in_file')]),
            (R2Star, R2star_neg_mask_reorient, [('neg_mask', 'in_file')]),
                                    
            #rename files and data sink
            (infosource, rename_infosource, [('subject_id', 'filename')]),
            #fieldmap
            (rename_infosource, rename_fieldmap, [('renamed','subject_id')]),            
            (fieldmap_reorient, rename_fieldmap, [('out_file', 'in_file')]),
            (rename_fieldmap, datasink, [('out_file', '@')]),
            #qsm
            (rename_infosource, rename_QSM, [('renamed','subject_id')]),                        
            (QSM_reorient, rename_QSM, [('out_file', 'in_file')]),
            (rename_QSM, datasink, [('out_file', '@.@qsm')]),
            #qsm brain mask
            (rename_infosource, rename_QSM_brain_mask, [('renamed','subject_id')]),                        
            (QSM_brain_mask_reorient, rename_QSM_brain_mask, [('out_file', 'in_file')]),
            (rename_QSM_brain_mask, datasink, [('out_file', '@.@qsm_brain')]),
            #qsm noisey voxels in fieldmap
            (rename_infosource, rename_QSM_noise_mask, [('renamed','subject_id')]),            
            (QSM_noise_mask_reorient, rename_QSM_noise_mask, [('out_file', 'in_file')]),
            (rename_QSM_noise_mask, datasink, [('out_file', '@.@qsm_noise')]),
            #r2star
            (rename_infosource, rename_R2star, [('renamed','subject_id')]),            
            (R2star_reorient, rename_R2star, [('out_file', 'in_file')]),
            (rename_R2star, datasink, [('out_file', '@.@r2star')]),   
            #r2star fit map
            (rename_infosource, rename_R2star_fit, [('renamed','subject_id')]),            
            (R2star_fit_reorient, rename_R2star_fit, [('out_file', 'in_file')]),
            (rename_R2star_fit, datasink, [('out_file', '@.@r2starfit')]), 
            #r2star negative values that were set to 0
            (rename_infosource, rename_R2star_neg_mask, [('renamed','subject_id')]),            
            (R2star_neg_mask_reorient, rename_R2star_neg_mask, [('out_file', 'in_file')]),
            (rename_R2star_neg_mask, datasink, [('out_file', '@.@r2starneg')]),
                        
            (infosource, datasink, [('subject_id', 'container')]),         
            ])
    return wf