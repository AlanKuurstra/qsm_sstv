import os
import unittest
import numpy as np
import nibabel as nib
import subprocess
from cfmm_nipype_interfaces.numerical_phantoms.qsmphantom import create_phantom
from cfmm_nipype_interfaces import MatlabRunMode
import tempfile


class TestQSMPipeline(unittest.TestCase):
    def __init__(self, *args, run_mode=MatlabRunMode.CONTAINER_MCR, **kwargs):
        self.run_mode = run_mode
        super(TestQSMPipeline, self).__init__(*args, **kwargs)

    def setUp(self):
        self.BidsDir_obj = tempfile.TemporaryDirectory()
        self.BidsDir = self.BidsDir_obj.name
        _, self.rois, self.susc_values_ppb, self.r2star_values, _, _ = create_phantom(bids_dir=self.BidsDir)

    def test_pipeline(self):
        matlab_executable = '/usr/local/matlab/R2016b/bin/matlab'
        mcr_location = '/storage/akuurstr/MATLAB_R2016b/v91/'

        run_executable = os.path.join(os.path.dirname(__file__), '..', 'run.py')
        outDir = tempfile.TemporaryDirectory().name
        parameters = [self.BidsDir, outDir, 'participant',
                      '--participant_label', 'FAKE',
                      '--SS_TV_lagrange_parameter', '0.01',
                      '--keep_unnecessary_outputs',
                      '--skip_fast',
                      '--mcr_location', mcr_location,
                      '--matlab_executable', matlab_executable,
                      '--run_mode', self.run_mode.name,
                      ]
        subprocess.run([run_executable, ] + parameters)

        # check to see if pipeline is accurate
        susc_file = os.path.join(outDir, 'qsm_sstv', 'sub-FAKE', 'anat', 'sub-FAKE_anat-QSM.nii.gz')
        susc = nib.load(susc_file).get_data()

        r2star_file = os.path.join(outDir, 'qsm_sstv', 'sub-FAKE', 'anat', 'sub-FAKE_anat-R2star.nii.gz')
        r2star = nib.load(r2star_file).get_data()

        for roi_index in range(1, self.rois.shape[-1]):
            susc_value_measured = susc[self.rois[..., roi_index].astype('bool')].mean()
            susc_value_true = self.susc_values_ppb[0] + self.susc_values_ppb[roi_index]
            self.assertLessEqual(np.abs((susc_value_measured - susc_value_true) / susc_value_true), 0.05,
                                 "ROI %s susceptibility failed. Measured: %s, Actual: %s" % (
                                     roi_index, susc_value_measured, susc_value_true))

            r2star_value_measured = r2star[self.rois[..., roi_index].astype('bool')].mean()
            r2star_value_true = self.r2star_values[0] + self.r2star_values[roi_index]
            self.assertLessEqual(np.abs((r2star_value_measured - r2star_value_true) / r2star_value_true), 0.05,
                                 "ROI %s r2 star failed. Measured: %s, Actual: %s" % (
                                     roi_index, r2star_value_measured, r2star_value_true))


if __name__ == '__main__':
    # make sure you run this file and not pycharm's "unittests for ..."  or "unittests in ..."
    suite = unittest.TestSuite()
    suite.addTest(TestQSMPipeline('test_pipeline', run_mode=MatlabRunMode.LOCAL_MCR))
    runner = unittest.TextTestRunner()
    runner.run(suite)
