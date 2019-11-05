#!/usr/bin/env python3
"""
@author: akuurstr

-in custom python nodes, using nibabel to open niftis and always copy the input
 affine matrix to the output nifti.
-in matlab nodes, using load_untouch_nii to load data without applying
 affine transformations. output nifti is saved using input's header
-user is expected to provide which dimension the nifti file stores the 
 B0 direction
"""
from nipype.interfaces import utility as niu
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.fsl as fsl
from bids import BIDSLayout
import os
from enum import Enum
from cfmm_nipype_interfaces.siemens_phase_preprocess import SiemensPhasePreprocess
from cfmm_nipype_interfaces.calc_r2star import CalcR2Star
from cfmm_nipype_interfaces.param_from_json import GetCFFromJson
from cfmm_nipype_interfaces.replace_slash import replace_slash
from cfmm_nipype_interfaces.calc_reliability_mask import CalculateReliabilityMask
from cfmm_nipype_interfaces.freq_from_phase import EstimateFrequencyFromWrappedPhase
from cfmm_nipype_interfaces.trim_mask import TrimMaskUsingReliability
from cfmm_nipype_interfaces.qsm.avg_and_snrmap import GetAvgAndSNRMap
from cfmm_nipype_interfaces.qsm.ss_tv import SS_TV

import cfmm_nipype_interfaces.qsm.matlab_scripts
matlab_scripts_path = cfmm_nipype_interfaces.qsm.matlab_scripts.__path__._path[0]




class RunMode(Enum):
    MATLAB = 1
    MCR = 2


class BrainExtractMethod(Enum):
    BET = 1
    BIDS = 2
    SINGLE_SUBJECT_FULL_PATH = 3


def create_pipeline_SS_TV(bids_dir,
                          work_dir,
                          out_dir,
                          subjects,
                          sessions,
                          mag_match_pattern,
                          phase_match_pattern,
                          mask_match_pattern,
                          keep_unnecessary_outputs,
                          FAST_bias_iters,
                          FAST_bias_lowpass,
                          FAST_num_classes,
                          skip_fast,
                          brain_extract_method,
                          BET_frac,
                          single_subject_custom_mask,
                          freq_weights__snr_window_sz,
                          truncate_echo,
                          SS_TV_lagrange_parameter,
                          B0_dir,
                          scnd_diff_reliability_thresh_noise,
                          trim_radius_sz,
                          scnd_diff_reliability_thresh_trim,
                          skip_qsm,
                          skip_r2star,
                          matlab_executable,
                          mcr_location,
                          run_mode):
    layout = BIDSLayout(bids_dir)

    ### CREATE PIPELINE OBJECT
    pipelineDir = work_dir
    wf = pe.Workflow(name="SS_TV")
    wf.base_dir = pipelineDir
    wf.config['execution']['remove_unnecessary_outputs'] = not keep_unnecessary_outputs

    ### GET MULTI-ECHO DATA
    # can we do this more elegantly?
    first_echo_files = []
    for subject in subjects:
        if layout.get_sessions(subject=subject) == []:
            if sessions == ['.*']:
                first_echo_files = first_echo_files + layout.get(subject=subject, modality='anat',
                                                                 extensions='.*part-phase.*echo-0*1.*.nii.*', )
            else:
                print(
                    "Warning: Session filter applied, but subject " + subject + " has no bids session information. This subject has been ignored.")
        else:
            for session in sessions:
                first_echo_files = first_echo_files + layout.get(subject=subject, session=session, modality='anat',
                                                                 extensions='.*part-phase.*echo-0*1.*.nii.*', )
    anat_folders = []
    for img in first_echo_files:
        full_dirname = os.path.dirname(img.filename)
        remove_base_dir = full_dirname.replace(bids_dir, '')
        remove_leading_slash = remove_base_dir.lstrip(os.sep)
        anat_folders.append(remove_leading_slash)

    anat_folders = list(set(anat_folders))
    anat_folders.sort()

    # IdentityInterface is useful for passing subject directory structure to datasink
    infosource = pe.Node(niu.IdentityInterface(fields=['subject_id']), name="infosource")
    infosource.iterables = ('subject_id', anat_folders)

    ### NODES AND PARAMETERS
    datasource = pe.Node(
        nio.DataGrabber(infields=['subject_id'],
                        outfields=['phase_images', 'mag_images', 'phase_jsons', 'mag_jsons', 'brain_mask']),
        name='datasource')
    datasource.inputs.field_template = dict(
        phase_images='%s/' + phase_match_pattern + '.nii*',
        phase_jsons='%s/' + phase_match_pattern + '.json',
        mag_images='%s/' + mag_match_pattern + '.nii*',
        mag_jsons='%s/' + mag_match_pattern + '.json',
        brain_mask='%s/' + mask_match_pattern + '.nii*',
    )
    # if brain_extract_method == BrainExtractMethod.BIDS:
    #    datasource.inputs.field_template['brain_mask'] = '%s/' + mask_match_pattern + '.nii*'
    datasource.inputs.sort_filelist = True
    datasource.inputs.template = "*"
    datasource.inputs.base_directory = bids_dir

    # this node must change depending on the scanner vendor
    susc_phase_preprocess = pe.Node(SiemensPhasePreprocess(), name='susc_phase_preprocess')

    avg_and_freq_estimate_weights = pe.Node(GetAvgAndSNRMap(), name='avg_and_freq_estimate_weights')
    avg_and_freq_estimate_weights.inputs.snr_window_sz = freq_weights__snr_window_sz
    avg_and_freq_estimate_weights.inputs.avg_out_filename = "avg.nii.gz"
    avg_and_freq_estimate_weights.inputs.snr_map_out_filename = "weights.nii.gz"

    wf.connect([
        (infosource, datasource, [('subject_id', 'subject_id')]),
        (datasource, avg_and_freq_estimate_weights, [('mag_images', 'mag')]),
        (datasource, susc_phase_preprocess, [('phase_images', 'infiles')])
    ])

    if brain_extract_method == BrainExtractMethod.BET:
        brain_extract = pe.Node(fsl.BET(), name='brain_extract_bet')
        brain_extract.inputs.frac = BET_frac
        brain_extract.inputs.mask = True
        brain_extract.inputs.robust = True

        if skip_fast:
            # connect avg directly to bet (skip FAST if image uniform enough for brain extraction)
            wf.connect([
                (avg_and_freq_estimate_weights, brain_extract, [('avg_out_filename', 'in_file')])
            ])
        else:
            # connect avg to nu correction, connect nu correction to bet
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
            nonuniformity_correct_fsl = pe.Node(fsl.FAST(), name='nonuniformity_correct_fsl')
            nonuniformity_correct_fsl.inputs.img_type = 2  # 1 for t1, 2 for t2
            nonuniformity_correct_fsl.inputs.bias_iters = FAST_bias_iters  # higher for larger nonuniformity
            nonuniformity_correct_fsl.inputs.bias_lowpass = FAST_bias_lowpass  # spm uses 30
            nonuniformity_correct_fsl.inputs.number_classes = FAST_num_classes  # spm uses 5
            nonuniformity_correct_fsl.inputs.output_biasfield = True
            nonuniformity_correct_fsl.inputs.output_biascorrected = True
            nonuniformity_correct_fsl.interface.estimated_memory_gb = 10

            wf.connect([
                # spm requires matlab
                # (avg_and_freq_estimate_weights, nonuniformityCorrect_spm, [('avgOutFilename', 'data')]),
                # (nonuniformityCorrect_spm, brain_extract, [('bias_corrected_image', 'in_file')]),
                (avg_and_freq_estimate_weights, nonuniformity_correct_fsl, [('avg_out_filename', 'in_files')]),
                (nonuniformity_correct_fsl, brain_extract, [('restored_image', 'in_file')])
            ])
    elif brain_extract_method == BrainExtractMethod.BIDS:
        brain_extract = pe.Node(
            nio.DataGrabber(infields=['subject_id'],
                            outfields=['mask_file']),
            name='bids_brain_mask')
        brain_extract.inputs.field_template = dict(
            mask_file='%s/' + mask_match_pattern + '.nii*',
        )
        brain_extract.inputs.sort_filelist = False
        brain_extract.inputs.template = "*"
        brain_extract.inputs.base_directory = bids_dir
        wf.connect([
            (infosource, brain_extract, [('subject_id', 'subject_id')]),
        ])

    elif brain_extract_method == BrainExtractMethod.SINGLE_SUBJECT_FULL_PATH:
        brain_extract = pe.Node(niu.IdentityInterface(fields=['mask_file']), name="fullpath_brain_mask")
        brain_extract.inputs.mask_file = single_subject_custom_mask

    freq_est = pe.Node(EstimateFrequencyFromWrappedPhase(), 'freq_est')
    freq_est.inputs.truncate_echo = truncate_echo
    freq_est.inputs.freq_filename = "freq_est.nii.gz"
    freq_est.interface.estimated_memory_gb = 4

    fieldmap_reorient = pe.Node(fsl.Reorient2Std(), name='fieldmap_reorient')

    datasink = pe.Node(nio.DataSink(), name="datasink")
    datasink.inputs.base_directory = out_dir + '/qsm_sstv/'
    datasink.inputs.parameterization = False

    rename_infosource = pe.Node(replace_slash, "rename_infosource")
    rename_fieldmap = pe.Node(niu.Rename(format_string="%(subject_id)s-fieldmap", keep_ext=True), "rename_fieldmap")

    wf.connect([
        (susc_phase_preprocess, freq_est, [('outfiles', 'phase')]),
        (datasource, freq_est, [('phase_jsons', 'json')]),
        (brain_extract, freq_est, [('mask_file', 'mask')]),
        (avg_and_freq_estimate_weights, freq_est, [('snr_map_out_filename', 'weight')]),
        (freq_est, fieldmap_reorient, [('freq_filename', 'in_file')]),
        # rename files and data sink
        (infosource, rename_infosource, [('subject_id', 'filename')]),
        # fieldmap
        (rename_infosource, rename_fieldmap, [('renamed', 'subject_id')]),
        (fieldmap_reorient, rename_fieldmap, [('out_file', 'in_file')]),
        (rename_fieldmap, datasink, [('out_file', '@')]),
        (infosource, datasink, [('subject_id', 'container')]),
    ])

    if not (skip_qsm and skip_r2star):
        trim_mask = pe.Node(TrimMaskUsingReliability(), name='trim_mask')
        trim_mask.inputs.erosion_sz = trim_radius_sz  # in mm
        trim_mask.inputs.threshold = scnd_diff_reliability_thresh_trim
        trim_mask.inputs.trimmed_mask_filename = "trim_mask.nii.gz"
        trim_mask.inputs.reliability_filename = "unreliableMap.nii.gz"
        trim_mask.interface.estimated_memory_gb = 25

        wf.connect([
            (freq_est, trim_mask, [('freq_filename', 'phase')]),
            (brain_extract, trim_mask, [('mask_file', 'mask')])
        ])

    if not skip_qsm:
        unreliable_fieldmap_voxels = pe.Node(CalculateReliabilityMask(), name='unreliable_fieldmap_voxels')
        unreliable_fieldmap_voxels.inputs.threshold = scnd_diff_reliability_thresh_noise
        unreliable_fieldmap_voxels.inputs.reliability_mask_filename = "unreliableMask.nii.gz"
        unreliable_fieldmap_voxels.inputs.reliability_filename = "unreliableMap.nii.gz"

        CF_value = pe.Node(GetCFFromJson, name='CFValue')

        matlab_execution = None
        if run_mode == RunMode.MATLAB:
            matlab_execution = matlab_executable + ' -nodisplay -nosplash -r "addpath ' + os.path.abspath(
                matlab_scripts_path) + ';SS_TV_script <params>;quit;"'
        elif run_mode == RunMode.MCR:
            matlab_execution = os.path.join(matlab_scripts_path,
                                            'run_SS_TV_script.sh') + ' ' + mcr_location + ' <params>'
        susceptibility = pe.Node(SS_TV(matlab_execution), name='susceptibility')
        susceptibility.inputs.quit_matlab = ''  # use this line when using mcr, comment when using matlab
        susceptibility.inputs.alpha = SS_TV_lagrange_parameter
        susceptibility.inputs.B0_dir = B0_dir
        susceptibility.inputs.susceptibility_filename = 'susceptibilityMap.nii.gz'
        susceptibility.interface.estimated_memory_gb = 10

        QSM_reorient = pe.Node(fsl.Reorient2Std(), name='QSM_reorient')
        QSM_brain_mask_reorient = pe.Node(fsl.Reorient2Std(), name='QSM_brain_mask_reorient')
        QSM_noise_mask_reorient = pe.Node(fsl.Reorient2Std(), name='QSM_noise_mask_reorient')

        rename_QSM = pe.Node(niu.Rename(format_string="%(subject_id)s-QSM", keep_ext=True), "rename_QSM")
        rename_QSM_brain_mask = pe.Node(niu.Rename(format_string="%(subject_id)s-QSM_brainMask", keep_ext=True),
                                        "rename_QSM_brain_mask")
        rename_QSM_noise_mask = pe.Node(niu.Rename(format_string="%(subject_id)s-QSM_noiseMask", keep_ext=True),
                                        "rename_QSM_noise_mask")
        wf.connect([
            (freq_est, unreliable_fieldmap_voxels, [('freq_filename', 'phase')]),
            (brain_extract, unreliable_fieldmap_voxels, [('mask_file', 'mask')]),
            (freq_est, susceptibility, [('freq_filename', 'freq_loc')]),
            (datasource, CF_value, [('mag_jsons', 'filename')]),
            (unreliable_fieldmap_voxels, susceptibility, [('reliability_mask_filename', 'reliability_mask_loc')]),
            (trim_mask, susceptibility, [('trimmed_mask_filename', 'mask_loc')]),
            (CF_value, susceptibility, [('CF_value', 'CF')]),

            (susceptibility, QSM_reorient, [('susceptibility_filename', 'in_file')]),
            (trim_mask, QSM_brain_mask_reorient, [('trimmed_mask_filename', 'in_file')]),
            (unreliable_fieldmap_voxels, QSM_noise_mask_reorient, [('reliability_mask_filename', 'in_file')]),

            # qsm
            (rename_infosource, rename_QSM, [('renamed', 'subject_id')]),
            (QSM_reorient, rename_QSM, [('out_file', 'in_file')]),
            (rename_QSM, datasink, [('out_file', '@.@qsm')]),
            # qsm brain mask
            (rename_infosource, rename_QSM_brain_mask, [('renamed', 'subject_id')]),
            (QSM_brain_mask_reorient, rename_QSM_brain_mask, [('out_file', 'in_file')]),
            (rename_QSM_brain_mask, datasink, [('out_file', '@.@qsm_brain')]),
            # qsm noisey voxels in fieldmap
            (rename_infosource, rename_QSM_noise_mask, [('renamed', 'subject_id')]),
            (QSM_noise_mask_reorient, rename_QSM_noise_mask, [('out_file', 'in_file')]),
            (rename_QSM_noise_mask, datasink, [('out_file', '@.@qsm_noise')]),
        ])

    if not skip_r2star:
        R2Star = pe.Node(CalcR2Star(), 'R2Star')
        R2Star.inputs.R2star = 'R2star.nii.gz'
        R2Star.inputs.neg_mask = 'negMask.nii.gz'
        R2Star.inputs.nan_mask = 'nanMask.nii.gz'
        # R2Star.interface.estimated_memory_gb = 5

        R2star_reorient = pe.Node(fsl.Reorient2Std(), name='R2star_reorient')
        R2star_fit_reorient = pe.Node(fsl.Reorient2Std(), name='R2star_fit_reorient')
        R2star_neg_mask_reorient = pe.Node(fsl.Reorient2Std(), name='R2star_neg_mask_reorient')

        rename_R2star = pe.Node(niu.Rename(format_string="%(subject_id)s-R2star", keep_ext=True), "rename_R2star")
        rename_R2star_fit = pe.Node(niu.Rename(format_string="%(subject_id)s-R2star_fit", keep_ext=True),
                                    "rename_R2star_fit")
        rename_R2star_neg_mask = pe.Node(niu.Rename(format_string="%(subject_id)s-R2star_negMask", keep_ext=True),
                                         "rename_R2star_neg_mask")

        wf.connect([
            (datasource, R2Star, [('mag_images', 'mag')]),
            (susc_phase_preprocess, R2Star, [('outfiles', 'phase')]),
            (freq_est, R2Star, [('freq_filename', 'freq_loc')]),
            (trim_mask, R2Star, [('trimmed_mask_filename', 'mask')]),
            (datasource, R2Star, [('mag_jsons', 'json')]),
            (R2Star, R2star_reorient, [('R2star', 'in_file')]),
            (R2Star, R2star_fit_reorient, [('R2star_fit', 'in_file')]),
            (R2Star, R2star_neg_mask_reorient, [('neg_mask', 'in_file')]),
            # r2star
            (rename_infosource, rename_R2star, [('renamed', 'subject_id')]),
            (R2star_reorient, rename_R2star, [('out_file', 'in_file')]),
            (rename_R2star, datasink, [('out_file', '@.@r2star')]),
            # r2star fit map
            (rename_infosource, rename_R2star_fit, [('renamed', 'subject_id')]),
            (R2star_fit_reorient, rename_R2star_fit, [('out_file', 'in_file')]),
            (rename_R2star_fit, datasink, [('out_file', '@.@r2starfit')]),
            # r2star negative values that were set to 0
            (rename_infosource, rename_R2star_neg_mask, [('renamed', 'subject_id')]),
            (R2star_neg_mask_reorient, rename_R2star_neg_mask, [('out_file', 'in_file')]),
            (rename_R2star_neg_mask, datasink, [('out_file', '@.@r2starneg')]),
        ])

    return wf
