#!/usr/bin/env python3
from pipeline_SS_TV_Siemens import create_pipeline_SS_TV, MatlabRunMode, BrainExtractMethod
import os
import argparse
from argparse import ArgumentParser
from nipype import config, logging

debugging = False

if __name__ == "__main__":
    # followed format from BIDS-Apps/nipypelines
    defstr = ' (default %(default)s)'
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('bids_dir', help='Data directory formatted according to BIDS standard.')
    parser.add_argument('output_dir', help='Directory where output files are to be stored.')
    parser.add_argument('analysis_level',
                        help='Level of the analysis that will be performed.',
                        choices=['participant'])
    parser.add_argument('--mag_match_pattern', dest="mag_match_pattern",
                        default='*part-mag_echo*',
                        help='Pattern used to match magnitude images and json files '
                             'in anat folder (leave extension out of pattern). The '
                             'pattern may contain simple shell-style wildcards a la '
                             'fnmatch. However, unlike fnmatch, filenames starting with '
                             'a dot are special cases that are not matched by \'*\' and '
                             '\'?\' patterns. Example usage: *part-mag_echo*', metavar=None)
    parser.add_argument('--phase_match_pattern', dest="phase_match_pattern",
                        default='*part-phase_echo*',
                        help='Pattern used to match phase images and json files '
                             'in anat folder (leave extension out of pattern). The '
                             'pattern may contain simple shell-style wildcards. '
                             'However, filenames starting with a dot are special cases '
                             'that are not matched by \'*\' and \'?\' patterns. '
                             'Example usage: *part-phase_echo*')
    parser.add_argument('--mask_match_pattern', dest="mask_match_pattern",
                        default='*brain_mask*',
                        help='Pattern used to match brain mask images '
                             'in BIDS anat folder (leave file extension out of pattern). '
                             'Must also set brain_extract_method parameter appropriately.The '
                             'pattern may contain simple shell-style wildcards. '
                             'However, filenames starting with a dot are special cases '
                             'that are not matched by \'*\' and \'?\' patterns. '
                             'Example usage: *brain_mask*')
    parser.add_argument('--participant_label',
                        help='The label(s) of the participant(s) that should be analyzed. The label '
                             'corresponds to sub-<participant_label> from the BIDS spec '
                             '(do not include prefix "sub-"). If this parameter is not '
                             'provided all subjects will be analyzed. Multiple '
                             'participants can be specified with a space separated list.',
                        default=['.*'],
                        nargs="+")
    parser.add_argument('--session_label', help='The label(s) of the session(s) that should be analyzed. The label '
                                                'corresponds to ses-<session_label> from the BIDS spec '
                                                '(do not include prefix "ses-"). If this parameter is not '
                                                'provided all sessions will be analyzed. Multiple '
                                                'sessions can be specified with a space separated list.',
                        default=['.*'],
                        nargs="+")
    parser.add_argument("--mcr_location", dest="mcr_location", default='/opt/mcr/v91',
                        help="Location of MATLAB Compiler Runtime.")
    parser.add_argument("-w", "--work_dir", dest="work_dir",
                        help="Work directory. Defaults to <output_dir>/scratch")
    parser.add_argument("-l", "--log_dir", dest="log_dir",
                        help="Nipype output log directory. Defaults to <output_dir>/log")
    parser.add_argument("-c", "--crash_dir", dest="crash_dir",
                        help="Nipype crash dump directory. Defaults to <output_dir>/crash_dump")
    parser.add_argument("-p", "--plugin", dest="plugin",
                        default='Linear', help="Nipype run plugin")
    parser.add_argument("--plugin_args", dest="plugin_args",
                        help="Nipype run plugin arguments")
    parser.add_argument("--keep_unnecessary_outputs", dest="keep_unnecessary_outputs",
                        action='store_true', default=False,
                        help="Keep all nipype node outputs, even if unused")
    parser.add_argument("--FAST_bias_iters", dest="FAST_bias_iters",
                        default=5, help="FAST bias iterations parameter ")
    parser.add_argument("--FAST_bias_lowpass", dest="FAST_bias_lowpass",
                        default=20, help="FAST bias lowpass parameter")
    parser.add_argument("--FAST_num_classes", dest="FAST_num_classes",
                        default=3, help="FAST number of classes parameter ")
    parser.add_argument("--skip_fast", dest="skip_fast",
                        action='store_true', default=False,
                        help="Skip FAST non-uniformity intensity correction. This saves time "
                             "but can affect brain extraction.")
    parser.add_argument("--brain_extract_method", choices=BrainExtractMethod.__members__, default='BET',
                        help="BET introduces an FSL BET node into the pipeline. BIDS finds a precomputed mask in the "
                             "BIDS directory using the mask_match_pattern parameter. SINGLE_SUBJECT_FULL_PATH finds "
                             "a precomputed mask using the single_subject_custom_mask parameter.")
    parser.add_argument("--BET_frac", dest="BET_frac",
                        default=0.2, help="BET frac parameter")
    parser.add_argument("--single_subject_custom_mask", dest="single_subject_custom_mask",
                        default=False,
                        help="Provide the full path to a custom brain mask in order to override BET node. "
                             "Must also set brain_extract_method appropriately. Only use this option when processing "
                             "a single subject.")
    parser.add_argument("--freq_weights_snr_window_sz", dest="freq_weights_snr__window_sz",
                        default=15.0, help="Size of local snr window, needed for weights in frequency fitting.")
    parser.add_argument("--truncate_echo", dest="truncate_echo",
                        default=0, help="How many echoes to include. Use -1 for all echoes")
    parser.add_argument("--SS_TV_lagrange_parameter", dest="SS_TV_lagrange_parameter",
                        default=0.35, help="Dipole inversion lagrange parameter for TV (higher value is more smooth)")
    parser.add_argument("--B0_dir", dest="B0_dir",
                        default=3, help='Value can be 1,2, or 3 and specifies which dimension '
                                        'in your data corresponds to the B0 direction')
    parser.add_argument("--scnd_diff_reliability_thresh_noise", dest="scnd_diff_reliability_thresh_noise",
                        default=1000000.0,
                        help='second difference reliability threshold - used to '
                             'mask fieldmap\'s noisey voxels (lower value masks more)')
    parser.add_argument("--trim_radius_sz", dest="trim_radius_sz",
                        default=15.0,
                        help='Radius (mm) of trim kernel searching boundary of brain mask to remove unreliable voxels. '
                             'Corrects liberal BET extractions.')
    parser.add_argument("--scnd_diff_reliability_thresh_trim", dest="scnd_diff_reliability_thresh_trim",
                        default=1000000.0,
                        help='second difference reliability threshold - used to trim fieldmap mask '
                             '(lower value trims more within trim kernel)')
    parser.add_argument("--skip_r2star", dest="skip_r2star",
                        action='store_true', default=False,
                        help='Skip R2* map calculation. Choose true if solely interested in the '
                             'unwrapped frequency map or susceptibility map.')
    parser.add_argument("--skip_qsm", dest="skip_qsm",
                        action='store_true', default=False,
                        help='Skip susceptibility map calculation. Choose true if solely interested '
                             'in the unwrapped frequency map or R2*.')
    parser.add_argument("--matlab_executable", dest="matlab_executable", default='/usr/local/matlab/R2016b/bin/matlab',
                        help=argparse.SUPPRESS)
    parser.add_argument("--run_mode", choices=MatlabRunMode.__members__, default=MatlabRunMode.CONTAINER_MCR.name,
                        help=argparse.SUPPRESS)

    if debugging:
        local_mcr = '/storage/akuurstr/MATLAB_R2016b/v91/'
        local_matlab = '/storage/akuurstr/MATLAB_R2016b/bin/matlab'
        BidsDir = os.path.abspath('test/test_data/bids')
        outDir = os.path.abspath('test/test_data/results')
        parameters = [BidsDir, outDir, 'participant',
                      '--participant_label', 'C001',
                      '--keep_unnecessary_outputs',
                      '--skip_fast',
                      '--mcr_location', local_mcr,
                      '--matlab_executable', local_matlab,
                      '--run_mode', 'MATLAB',
                      ]
        args = parser.parse_args(parameters)
    else:
        args = parser.parse_args()

    bids_dir = args.bids_dir
    out_dir = args.output_dir

    mag_match_pattern = args.mag_match_pattern
    phase_match_pattern = args.phase_match_pattern
    mask_match_pattern = args.mask_match_pattern
    subjects = args.participant_label
    sessions = args.session_label

    if args.work_dir:
        work_dir = os.path.abspath(args.work_dir)
    else:
        work_dir = os.path.join(out_dir, 'scratch')
    if args.log_dir:
        log_dir = os.path.abspath(args.log_dir)
    else:
        tmp = "log-" + "_".join(subjects) + '-' + "_".join(sessions)
        tmp = tmp.replace(".*", "all").replace("*", "star")
        log_dir = os.path.join(out_dir, 'logs', tmp)
    if args.crash_dir:
        crash_dir = os.path.abspath(args.crash_dir)
    else:
        crash_dir = os.path.join(out_dir, 'crash_dump')
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    config.update_config({'logging': {
        'log_directory': log_dir,
        'log_to_file': True,
    },
        'execution': {
            'crashdump_dir': crash_dir,
            'crashfile_format': 'txt',
        }})
    logging.update_logging(config)

    plugin = args.plugin
    plugin_args = args.plugin_args
    keep_unnecessary_outputs = args.keep_unnecessary_outputs

    FAST_bias_iters = int(args.FAST_bias_iters)
    FAST_bias_lowpass = int(args.FAST_bias_lowpass)
    FAST_num_classes = int(args.FAST_num_classes)
    skip_fast = args.skip_fast
    brain_extract_method = BrainExtractMethod[args.brain_extract_method]
    BET_frac = float(args.BET_frac)
    single_subject_custom_mask = args.single_subject_custom_mask
    freq_weights__snr_window_sz = float(args.freq_weights_snr__window_sz)
    truncate_echo = int(args.truncate_echo)
    SS_TV_lagrange_parameter = float(args.SS_TV_lagrange_parameter)
    B0_dir = int(args.B0_dir)
    scnd_diff_reliability_thresh_noise = float(args.scnd_diff_reliability_thresh_noise)
    trim_radius_sz = float(args.trim_radius_sz)
    scnd_diff_reliability_thresh_trim = float(args.scnd_diff_reliability_thresh_trim)
    skip_qsm = args.skip_qsm
    skip_r2star = args.skip_r2star
    matlab_executable = args.matlab_executable
    mcr_location = args.mcr_location
    run_mode = MatlabRunMode[args.run_mode]

    wf_SS_TV = create_pipeline_SS_TV(bids_dir,
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
                                     run_mode
                                     )
    if args.plugin_args:
        execGraph_SS_TV = wf_SS_TV.run(args.plugin, plugin_args=eval(args.plugin_args))
    else:
        execGraph_SS_TV = wf_SS_TV.run(args.plugin)
