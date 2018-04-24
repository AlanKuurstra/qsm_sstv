#!/usr/bin/env python
from pipeline_SS_TV_Siemens import create_pipeline_SS_TV
import os

if __name__ == "__main__":
    
    #followed format from BIDS-Apps/nipypelines
    from argparse import ArgumentParser, RawTextHelpFormatter
    from nipype import config, logging
    defstr = ' (default %(default)s)'
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument('bids_dir', help='The directory with the input dataset '
                        'formatted according to the BIDS standard.')
    parser.add_argument('output_dir', help='The directory where the output files '
                        'should be stored. If you are running group level analysis '
                        'this folder should be prepopulated with the results of the'
                        'participant level analysis.')
    parser.add_argument('analysis_level', help='Level of the analysis that will be performed. '
                        'Multiple participant level analyses can be run independently '
                        '(in parallel) using the same output_dir.',
                        choices=['participant'], nargs='?', default='participant')
    parser.add_argument('--mag_match_pattern', dest="mag_match_pattern",
                        default='*part-mag*',
                        help='Pattern used to match magnitude images and json files '
                        'in anat folder (leave extension out of pattern). The '
                        'pattern may contain simple shell-style wildcards a la '
                        'fnmatch. However, unlike fnmatch, filenames starting with '
                        'a dot are special cases that are not matched by \'*\' and '
                        '\'?\' patterns. Example usage: *part-mag*')  
    parser.add_argument('--phase_match_pattern', dest="phase_match_pattern",
                        default='*part-phase*',
                        help='Pattern used to match phase images and json files '
                        'in anat folder (leave extension out of pattern). The '
                        'pattern may contain simple shell-style wildcards. '
                        'However, filenames starting with a dot are special cases '
                        'that are not matched by \'*\' and \'?\' patterns. '
                        'Example usage: *part-phase*')  
    parser.add_argument('--participant_label', help='The label(s) of the participant(s) that should be analyzed. The label '
                        'corresponds to sub-<participant_label> from the BIDS spec '
                        '(so it does not include "sub-"). If this parameter is not '
                        'provided all subjects should be analyzed. Multiple '
                        'participants can be specified with a space separated list.',
                        default=['.*'],
                        nargs="+")    
    parser.add_argument('--session_label', help='The label(s) of the session(s) that should be analyzed. The label '
                        'corresponds to ses-<session_label> from the BIDS spec '
                        '(so it does not include "ses-"). If this parameter is not '
                        'provided all sessions should be analyzed. Multiple '
                        'sessions can be specified with a space separated list.',
                        default=['.*'],
                        nargs="+")    
    parser.add_argument("-w", "--work_dir", dest="work_dir",
                        help="Work directory. Defaults to <output_dir>/scratch")
    parser.add_argument("-l", "--log_dir", dest="log_dir",
                        help="Nipype output log directory. Defaults to <output_dir>/log")
    parser.add_argument("-c", "--crash_dir", dest="crash_dir",
                        help="Nipype crash dump directory. Defaults to <output_dir>/crash_dump")
    parser.add_argument("-p", "--plugin", dest="plugin",
                        default='Linear',
                        help="Plugin to use")
    parser.add_argument("--plugin_args", dest="plugin_args",
                        help="Plugin arguments")
    parser.add_argument("--keep_unnecessary_outputs", dest="keep_unnecessary_outputs",
                        action='store_true',default=False,
                        help="keep all nipype node outputs, even if unused")
    parser.add_argument("--FAST_bias_iters", dest="FAST_bias_iters",
                        default=5,
                        help="FAST bias iterations parameter ")    
    parser.add_argument("--FAST_bias_lowpass", dest="FAST_bias_lowpass",
                        default=20,
                        help="FAST bias lowpass parameter")
    parser.add_argument("--FAST_num_classes", dest="FAST_num_classes",
                        default=3,
                        help="FAST number of classes parameter ")    
    parser.add_argument("--BET_frac", dest="BET_frac",
                        default=0.2,
                        help="BET frac parameter")
    parser.add_argument("--freq_weights_snr_window_sz", dest="freq_weights_snr__window_sz",
                        default=15.0,
                        help="Size of local snr window, needed for weights in frequency fitting.")
    parser.add_argument("--truncate_echo", dest="truncate_echo",
                        default=0,
                        help="How many echoes to include. Use -1 for all echoes")
    parser.add_argument("--SS_TV_lagrange_parameter", dest="SS_TV_lagrange_parameter",
                        default=0.35,
                        help="Dipole inversion lagrange parameter for TV")
    parser.add_argument("--B0_dir", dest="B0_dir",
                        default=3,
                        help="Value can be 1,2, or 3 and specifies which dimension '
                        'in your data corresponds to the B0 direction")
    parser.add_argument("--scnd_diff_reliability_thresh_trim", dest="scnd_diff_reliability_thresh_trim",
                        default=1000000.0,
                        help="second difference reliability threshold - used to trim fieldmap")
    parser.add_argument("--scnd_diff_reliability_thresh_noise", dest="scnd_diff_reliability_thresh_noise",
                        default=1000000.0,
                        help="second difference reliability threshold - used to '
                        'mask fieldmap's noisey voxels")
    
    args = parser.parse_args()
    
    #for testing
    #BidsDir='/workspace/akuurstr/ali_khan'
    #outDir='/workspace/akuurstr/ali_khan_results'     
    #args = parser.parse_args([BidsDir, outDir,'--participant_label','C011','--mag_match_pattern','*part-mag_echo*',
    #      '--phase_match_pattern','*part-phase_echo*','--SS_TV_lagrange_parameter','0.4','--truncate_echo','3',
    #      '--keep_unnecessary_outputs','--B0_dir','1'])    
 
    bids_dir=args.bids_dir
    out_dir=args.output_dir
    
    mag_match_pattern=args.mag_match_pattern
    phase_match_pattern=args.phase_match_pattern  
    subjects=args.participant_label
    sessions=args.session_label
        
    if args.work_dir:
        work_dir = os.path.abspath(args.work_dir)
    else:
        work_dir = os.path.join(out_dir, 'scratch')
    if args.log_dir:
        log_dir = os.path.abspath(args.log_dir)
    else:
        tmp="log-"+"_".join(subjects)+'-'+"_".join(sessions)
        tmp=tmp.replace(".*","all").replace("*","star")
        log_dir = os.path.join(out_dir, 'logs',tmp)
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
    
    plugin=args.plugin
    plugin_args=args.plugin_args
    keep_unnecessary_outputs=args.keep_unnecessary_outputs
        
    FAST_bias_iters=int(args.FAST_bias_iters)
    FAST_bias_lowpass=int(args.FAST_bias_lowpass)
    FAST_num_classes=int(args.FAST_num_classes)
    BET_frac=float(args.BET_frac)
    freq_weights__snr_window_sz=float(args.freq_weights_snr__window_sz)
    truncate_echo=int(args.truncate_echo)
    SS_TV_lagrange_parameter=float(args.SS_TV_lagrange_parameter)
    B0_dir=int(args.B0_dir)    
    scnd_diff_reliability_thresh_trim=float(args.scnd_diff_reliability_thresh_trim)
    scnd_diff_reliability_thresh_noise=float(args.scnd_diff_reliability_thresh_noise)
    
    wf_SS_TV=create_pipeline_SS_TV(bids_dir,
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
                                    scnd_diff_reliability_thresh_noise
                                   )
    if args.plugin_args:
        execGraph_SS_TV=wf_SS_TV.run(args.plugin, plugin_args=eval(args.plugin_args))
    else:
        execGraph_SS_TV=wf_SS_TV.run(args.plugin)    
