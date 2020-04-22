import os
import shutil
import subprocess
import unittest

import pandas as pd
from testfixtures import TempDirectory

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
PIPELINE_BASE_DIR = os.path.abspath(os.path.join(TEST_DIR, '..'))

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
        with TempDirectory() as temp_dir:
            os.chdir(temp_dir.path)

            example_data_dir = os.path.join(PIPELINE_BASE_DIR, 'data')
            results_dir = os.path.join(temp_dir.path, 'results')
            atac_general_config = os.path.join(PIPELINE_BASE_DIR, 'config', 'ATAC_general.yaml')

            ########################################
            # Create the configfile
            config_creator_script = os.path.join(PIPELINE_BASE_DIR, 'scripts', 'ngs_rawdata_config_creator.py')
            general_input = atac_general_config
            per_lib_input = os.path.join(example_data_dir, 'atac_test_data', 'atac_test_samplesheet.csv')
            results_dir = results_dir
            results_temp_dir = os.path.join(results_dir, 'tmp')
            outputfile_config = os.path.join(temp_dir.path, 'config_test.yaml')

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
            dryrun_command_fmt = 'snakemake --snakefile {snakefile} --configfile {configfile} -np'

            dryrun_command = dryrun_command_fmt.format(
                snakefile = os.path.join(PIPELINE_BASE_DIR, 'Snakefile_ATACseq'),
                configfile = os.path.join(temp_dir.path, 'config_test.yaml')
            )

            dryrun_return_code = subprocess.call(dryrun_command, shell=True)

            self.assertEqual(0, dryrun_return_code)

    # def test_endtoend_passes(self):
    #     with TempDirectory() as temp_dir:
    #         os.chdir(temp_dir.path)
    #
    #         example_data_dir = os.path.join(PIPELINE_BASE_DIR, 'data')
    #         results_dir = os.path.join(temp_dir.path, 'results')
    #         atac_general_config = os.path.join(PIPELINE_BASE_DIR, 'config', 'ATAC_general.yaml')
    #         cluster_config = os.path.join(PIPELINE_BASE_DIR, 'config', 'cluster_config.yaml')
    #
    #         ########################################
    #         # Create the configfile
    #         cfg_create_command_fmt = '{config_creator_script} --general_input {general_input} --per_lib_input {per_lib_input} --results_dir {results_dir} --temp_dir {temp_dir} > {outputfile_config}'
    #         cfg_create_command = cfg_create_command_fmt.format(
    #             config_creator_script = os.path.join(PIPELINE_BASE_DIR, 'scripts', 'ngs_rawdata_config_creator.py'),
    #             general_input = atac_general_config,
    #             per_lib_input = os.path.join(example_data_dir, 'atac_test_data', 'atac_test_samplesheet.csv'),
    #             results_dir = results_dir,
    #             temp_dir = os.path.join(results_dir, 'tmp'),
    #             outputfile_config = os.path.join(results_dir, 'config_test.yaml')
    #         )
    #         try:
    #             #return code only available from subprocess.check_output if non-zero (raises CalledProcessError)
    #             cfg_create_return_code = 0
    #             actual_output = subprocess.check_output(cfg_create_command, shell=True)
    #             # actual_output = actual_output.decode("utf-8")
    #             # lines = actual_output.split('\n')
    #             # empty_last_line = lines.pop()
    #             # dryrun_line = lines.pop()
    #         except subprocess.CalledProcessError as e:
    #             cfg_create_return_code = e.returncode
    #             output = e.output.decode()
    #             print(output)
    #
    #         self.assertEqual(0, cfg_create_return_code)
    #
    #         ########################################
    #         # End-to-end test
    #         endtoend_command_fmt = 'snakemake --snakefile {snakefile} --configfile {configfile} --latency-wait 60 --jobs 144 --cluster-config {cluster_config} --use-conda --cluster "sbatch --job-name={{cluster.name}} --account={{cluster.account}} --partition={{cluster.partition}} --nodes={{cluster.nodes}} --ntasks-per-node={{cluster.ntask}} --mem={{cluster.memory}} --time={{cluster.time}} --output=logs/%x-%j.out"'
    #
    #         endtoend_command = endtoend_command_fmt.format(
    #             snakefile = os.path.join(PIPELINE_BASE_DIR, 'Snakefile_ATACseq'),
    #             configfile = os.path.join(results_dir, 'config_test.yaml'),
    #             cluster_config = cluster_config
    #         )
    #
    #         return_code = subprocess.call(endtoend_command, shell=True)
    #
    #         self.assertEqual(1, 2)
