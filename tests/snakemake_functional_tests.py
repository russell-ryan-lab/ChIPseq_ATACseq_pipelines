import os
import shutil
import subprocess
import tempfile
import unittest

import pandas as pd

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
PIPELINE_BASE_DIR = os.path.abspath(os.path.join(TEST_DIR, '..'))
TEST_TMP_DIR = '/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/tmp'

class SnakemakeFunctionalTests(unittest.TestCase):
    def setUp(self):
        self.original_wd = os.getcwd()
        # self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqualEnough)

    def tearDown(self):
        os.chdir(self.original_wd)

    def assertDataframeEqualEnough(self, a, b, msg):
        #https://stackoverflow.com/a/54344148
        # allows us to use pandas testing.assert_frame_equal within unittest.testcase
        # Note using check_exact=False
        try:
            pd.testing.assert_frame_equal(a, b, check_exact=False)
        except AssertionError as e:
            raise self.failureException(msg) from e


    def test_dryrun_passes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            os.chdir(temp_dir)

            example_data_dir = os.path.join(PIPELINE_BASE_DIR, 'data')
            results_dir = os.path.join(temp_dir, 'results')
            results_temp_dir = os.path.join(results_dir, 'tmp')
            atac_general_config = os.path.join(PIPELINE_BASE_DIR, 'config', 'ATAC_general.yaml')

            ########################################
            # Create the configfile
            config_creator_script = os.path.join(PIPELINE_BASE_DIR, 'scripts', 'config_creator.py')
            general_input = atac_general_config
            per_lib_input = os.path.join(example_data_dir, 'atac_test_data', 'atac_test_samplesheet.csv')
            outputfile_config = os.path.join(temp_dir, 'config_test.yaml')

            with open(outputfile_config, 'w') as outfile_handle:
                cfg_create_return_code = subprocess.call([
                    config_creator_script,
                    '--general_input', general_input,
                    '--per_lib_input', per_lib_input,
                    '--results_dir', results_dir,
                    '--temp_dir', results_temp_dir
                ], stdout=outfile_handle)

                self.assertEqual(0, cfg_create_return_code)

            ########################################
            # Dryrun test
            snakefile = os.path.join(PIPELINE_BASE_DIR, 'ATACseq.smk')
            configfile = os.path.join(temp_dir, 'config_test.yaml')

            dryrun_return_code = subprocess.call([
                'snakemake',
                '--snakefile', snakefile,
                '--configfile', configfile,
                '-np'
            ])

            self.assertEqual(0, dryrun_return_code)

    def test_endtoend_passes(self):
        # Using specific temporary directory prefix so that worker nodes on GreatLakes cluster
        # can all access the same temp space. If this is not done, TemporaryDirectory creates them
        # in /tmp, which is not shared between nodes
        with tempfile.TemporaryDirectory(prefix = TEST_TMP_DIR) as temp_dir:
            os.chdir(temp_dir)

            example_data_dir = os.path.join(PIPELINE_BASE_DIR, 'data')
            results_dir = os.path.join(temp_dir, 'results')
            results_temp_dir = os.path.join(results_dir, 'tmp')
            results_log_dir = os.path.join(results_dir, 'logs')
            atac_general_config = os.path.join(PIPELINE_BASE_DIR, 'config', 'ATAC_general.yaml')

            # Create /logs and /tmp folder in results_dir
            os.mkdir(results_dir)
            os.mkdir(results_temp_dir)
            os.mkdir(results_log_dir)

            ########################################
            # Create the configfile
            config_creator_script = os.path.join(PIPELINE_BASE_DIR, 'scripts', 'config_creator.py')
            general_input = atac_general_config
            per_lib_input = os.path.join(example_data_dir, 'atac_test_data', 'atac_test_samplesheet.csv')
            outputfile_config = os.path.join(temp_dir, 'config_test.yaml')

            with open(outputfile_config, 'w') as outfile_handle:
                cfg_create_return_code = subprocess.call([
                    config_creator_script,
                    '--general_input', general_input,
                    '--per_lib_input', per_lib_input,
                    '--results_dir', results_dir,
                    '--temp_dir', results_temp_dir
                ], stdout=outfile_handle)

                self.assertEqual(0, cfg_create_return_code)

            ########################################
            # End-to-end test
            snakefile = os.path.join(PIPELINE_BASE_DIR, 'ATACseq.smk')
            configfile = os.path.join(temp_dir, 'config_test.yaml')
            slurm_status_script = os.path.join(PIPELINE_BASE_DIR, 'scripts', 'slurm_status.py')
            cluster_config = os.path.join(PIPELINE_BASE_DIR, 'config', 'cluster_config.yaml')

            endtoend_return_code = subprocess.call([
                'snakemake',
                '--snakefile', snakefile,
                '--configfile', configfile,
                '--latency-wait', '60',
                '--jobs', '144',
                '--cluster-config', cluster_config,
                '--use-conda',
                '--conda-prefix', '/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/atac_test_run2/.snakemake/conda',
                '--cluster-status', slurm_status_script,
                '--cluster',
                'sbatch --parsable --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --output=logs/%x-%j.out'
            ])

            self.assertEqual(0, endtoend_return_code)
