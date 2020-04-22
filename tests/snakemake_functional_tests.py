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
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqualEnough)

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
            tmp_data_dir = os.path.join(temp_dir.path, 'results')
            #Copy example data here
            atac_general_config = os.path.join(PIPELINE_BASE_DIR, 'config', 'ATAC_general.yaml')
            shutil.copytree(example_data_dir, tmp_data_dir)
            shutil.copyfile(atac_general_config, os.path.join(tmp_data_dir, 'ATAC_general.yaml'))
            #Create the configfile
            cfg_create_command_fmt = '{config_creator_script} --general_input {general_input} --per_lib_input {per_lib_input} --results_dir {results_dir} --temp_dir {temp_dir} > {outputfile_config}'
            cfg_create_command = cfg_create_command_fmt.format(
                config_creator_script = os.path.join(PIPELINE_BASE_DIR, 'scripts', 'ngs_rawdata_config_creator.py'),
                general_input = atac_general_config,
                per_lib_input = os.path.join(example_data_dir, 'atac_test_data', 'atac_test_samplesheet.csv'),
                results_dir = tmp_data_dir,
                temp_dir = os.path.join(tmp_data_dir, 'tmp'),
                outputfile_config = os.path.join(tmp_data_dir, 'config_test.yaml')
            )
            try:
                #return code only available from subprocess.check_output if non-zero (raises CalledProcessError)
                cfg_create_return_code = 0
                actual_output = subprocess.check_output(cfg_create_command, shell=True)
                # actual_output = actual_output.decode("utf-8")
                # lines = actual_output.split('\n')
                # empty_last_line = lines.pop()
                # dryrun_line = lines.pop()
            except subprocess.CalledProcessError as e:
                cfg_create_return_code = e.returncode
                output = e.output.decode()
                print(output)

            self.assertEqual(0, cfg_create_return_code)

            dryrun_command_fmt = 'snakemake --snakefile {snakefile} --configfile {configfile} -np'

            dryrun_command = dryrun_command_fmt.format(
                snakefile = os.path.join(PIPELINE_BASE_DIR, 'Snakefile_ATACseq'),
                configfile = os.path.join(tmp_data_dir, 'config_test.yaml')
            )
            try:
                #return code only available from subprocess.check_output if non-zero (raises CalledProcessError)
                dryrun_return_code = 0
                actual_output = subprocess.check_output(dryrun_command, shell=True)
                actual_output = actual_output.decode("utf-8")
                lines = actual_output.split('\n')
                empty_last_line = lines.pop()
                dryrun_line = lines.pop()
            except subprocess.CalledProcessError as e:
                dryrun_return_code = e.returncode
                output = e.output.decode()
                print(output)

            self.assertEqual(0, dryrun_return_code)


# --latency-wait 60 --jobs 144 --cluster-config {cluster_config} --use-conda --cluster "sbatch --job-name={cluster.name} --account={{cluster.account}} --partition={{cluster.partition}} --nodes={{cluster.nodes}} --ntasks-per-node={{cluster.ntask}} --mem={{cluster.memory}} --time={{cluster.time}} --output=logs/%x-%j.out"



        #     #Copy samplesheet here
        #     shutil.copyfile(os.path.join(WATERMELON_BASE_DIR, 'config', 'example_samplesheet.csv'), 'samplesheet.csv')
        #
        #     config_replacement_vals = {
        #         'samplesheet': os.path.join(temp_dir.path, 'samplesheet.csv'),
        #         'dirs': {'input' : os.path.join(tmp_data_dir, 'sim_reads_human_chr22')},
        #         'references': {
        #             'fasta' : os.path.join(tmp_data_dir, 'Homo_sapiens.GRCh38.dna_sm.chr22.fa.gz'), # These just need to exist for dry-run
        #             'gtf' : os.path.join(tmp_data_dir, 'Homo_sapiens.GRCh38.98.chr22.gtf.gz'), # For real-deal, unzip them first
        #             'annotation_tsv' : os.path.join(tmp_data_dir, 'Homo_sapiens.GRCh38.98_annotation.tsv.gz') # then point to those non-gz files
        #         },
        #         'email' : {
        #             'to': 'nobody'
        #         },
        #     }
        #
        #     #Create modified config in this temp dir, using example config and replacing values as needed
        #     configfile_path = os.path.join(temp_dir.path, 'testcase_config.yaml')
        #     testing_utils.create_modified_config(EXAMPLE_CONFIGFILE_PATH, configfile_path, config_replacement_vals)
        #
        #     command_fmt = 'snakemake --snakefile {} --configfile {} -n --config skip_validation=True {}'
        #     command = command_fmt.format(SNAKEFILE_PATH, configfile_path, REDIRECT_OUTPUT)
        #     print(command)
        #     try:
        #         #return code only available from subprocess.check_output if non-zero (raises CalledProcessError)
        #         return_code = 0
        #         actual_output = subprocess.check_output(command, shell=True)
        #         actual_output = actual_output.decode("utf-8")
        #         lines = actual_output.split('\n')
        #         empty_last_line = lines.pop()
        #         dryrun_line = lines.pop()
        #     except subprocess.CalledProcessError as e:
        #         return_code = e.returncode
        #         output = e.output.decode()
        #         print(output)
        #
        # self.assertEqual(0, return_code)
        # self.assertRegex(dryrun_line, 'This was a dry-run')
