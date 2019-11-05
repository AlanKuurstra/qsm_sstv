import os
import unittest
import numpy as np
import nibabel as nib
import hashlib
import pickle as pkl
import subprocess
import shutil

# ======================================================================================================
# only set this to true when you're sure the pipeline is outputting what you want and you need to update
# the hash values.
save_new_md5_dict = False
# ======================================================================================================

md5_dict_filename = os.path.join(os.path.dirname(__file__), 'qsm_output_md5_values.pkl')
if save_new_md5_dict:
    md5_dict = {}
else:
    with open(md5_dict_filename, 'rb') as f:
        md5_dict = pkl.load(f)


def hash_scientific_rounded(numpy_array):
    # need to truncate the float mantissa due to floating point errors differing between local system that builds
    # the checksum dict and the docker image that is running continuous integration unit testing
    sig = hashlib.md5()
    for value in numpy_array.ravel():
        scientific_format = np.format_float_scientific(value, precision=8)
        sig.update(scientific_format.encode())
    return sig.digest()


def find_nifti_images(node_outvalue_list, nifti_location_list, strip_string=None, compute_md5=False):
    for list_element in node_outvalue_list:
        if type(list_element) == list:
            find_nifti_images(list_element, nifti_location_list, compute_md5)
        elif type(list_element) == str:
            img_loc = list_element
            if os.path.exists(img_loc):
                try:
                    nifti_obj = nib.load(img_loc)
                    nifti_location_list.append(img_loc)
                    if compute_md5:
                        # need to strip the output_dir from the key because it will likely be different when the
                        # dictionary is generated on a local machine and when it is checked in the gitlab-runner
                        # docker environment
                        if strip_string:
                            img_loc_stripped = img_loc.replace(strip_string, '')
                        else:
                            img_loc_stripped = img_loc
                        md5_dict[img_loc_stripped] = hash_scientific_rounded(nifti_obj.get_data())
                except:
                    print(img_loc, "is not a nifti!")


class TestQSMPipeline(unittest.TestCase):
    def setUp(self):
        self.local_mcr = '/storage/akuurstr/MATLAB_R2016b/v91/'
        self.BidsDir = os.path.abspath('test_data/bids')
        self.outDir = os.path.abspath('test_data/results')

    def test_calculation(self):
        run_executable = os.path.abspath('../run.py')
        parameters = [self.BidsDir, self.outDir, 'participant',
                      '--participant_label', 'C001',
                      '--keep_unnecessary_outputs',
                      '--skip_fast',
                      '--mcr_location', self.local_mcr,
                      ]
        subprocess.run([run_executable, ] + parameters)

        nifti_output_files = []
        all_output_files = []
        for root, dirs, files in os.walk(os.path.join(self.outDir, 'qsm_sstv')):
            for file in files:
                all_output_files.append(file)
        find_nifti_images(all_output_files, nifti_output_files, self.outDir, save_new_md5_dict)

        for img_loc in nifti_output_files:
            img_loc_stripped = img_loc.replace(self.outDir, '')
            self.assertEqual(md5_dict[img_loc_stripped],
                             hash_scientific_rounded(nib.load(img_loc).get_data()),
                             img_loc_stripped + " does not match hashed value.")

        if save_new_md5_dict:
            with open(md5_dict_filename, 'wb') as f:
                pkl.dump(md5_dict, f)

    def tearDown(self):
        shutil.rmtree(self.outDir)


if __name__ == '__main__':
    unittest.main()
