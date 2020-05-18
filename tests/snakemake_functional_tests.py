import os
import shutil
import subprocess
import tempfile
import unittest

import pandas as pd

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
PIPELINE_BASE_DIR = os.path.abspath(os.path.join(TEST_DIR, '..'))
TEST_TMP_DIR = os.path.join(PIPELINE_BASE_DIR, '')

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


    def test_atac_dryrun_passes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            os.chdir(temp_dir)

            results_dir = os.path.join(temp_dir, 'results')
            results_temp_dir = os.path.join(results_dir, 'tmp')
            atac_general_config = os.path.join(PIPELINE_BASE_DIR, 'config', 'ATAC_general.yaml')

            ########################################
            # Create the configfile
            config_creator_script = os.path.join(PIPELINE_BASE_DIR, 'scripts', 'config_creator.py')
            general_input = atac_general_config

            test_samplesheet = '''lib,sample,genome,basepath
1,s1_control,mm10,{}/data/atac_test_data/Sample_1
2,s1_treatment,mm10,{}/data/atac_test_data/Sample_2
3,s2_treatment,mm10,{}/data/atac_test_data/Sample_3
4,s3_treatment,mm10,{}/data/atac_test_data/Sample_4'''
            test_samplesheet = test_samplesheet.format(PIPELINE_BASE_DIR, PIPELINE_BASE_DIR, PIPELINE_BASE_DIR, PIPELINE_BASE_DIR)
            test_sheet_filename = os.path.join(temp_dir, 'samplesheet_test.csv')
            with open(test_sheet_filename, 'w') as test_sheet_filehandle:
                test_sheet_filehandle.write(test_samplesheet)

            outputfile_config = os.path.join(temp_dir, 'config_test.yaml')

            with open(outputfile_config, 'w') as outfile_handle:
                cfg_create_return_code = subprocess.call([
                    config_creator_script,
                    '--general_input', general_input,
                    '--per_lib_input', test_sheet_filename,
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

    def test_chip_se_dryrun_passes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            os.chdir(temp_dir)

            results_dir = os.path.join(temp_dir, 'results')
            results_temp_dir = os.path.join(results_dir, 'tmp')
            chip_general_config = os.path.join(PIPELINE_BASE_DIR, 'config', 'ChIP_histone_general_se.yaml')

            ########################################
            # Create the configfile
            config_creator_script = os.path.join(PIPELINE_BASE_DIR, 'scripts', 'config_creator.py')
            general_input = chip_general_config

            test_samplesheet = '''lib,sample,genome,input,homer_fmg_genome,basepath
SRR6214070,GM12878_NKRF,hg19,GM12878_Input,hg19r,{}/data/sra_chip_test_data/SAMN06121956/
SRR5112307,GM12878_Input,hg19,,,{}/data/sra_chip_test_data/SAMN06122176/''' # FIXME: make these samples match new chip test data when we have it
            test_samplesheet = test_samplesheet.format(PIPELINE_BASE_DIR, PIPELINE_BASE_DIR)
            test_sheet_filename = os.path.join(temp_dir, 'samplesheet_test.csv')
            with open(test_sheet_filename, 'w') as test_sheet_filehandle:
                test_sheet_filehandle.write(test_samplesheet)

            outputfile_config = os.path.join(temp_dir, 'config_test.yaml')

            with open(outputfile_config, 'w') as outfile_handle:
                cfg_create_return_code = subprocess.call([
                    config_creator_script,
                    '--general_input', general_input,
                    '--per_lib_input', test_sheet_filename,
                    '--results_dir', results_dir,
                    '--temp_dir', results_temp_dir
                ], stdout=outfile_handle)

                self.assertEqual(0, cfg_create_return_code)

            ########################################
            # Dryrun test
            snakefile = os.path.join(PIPELINE_BASE_DIR, 'ChIPseq_se.smk')
            configfile = os.path.join(temp_dir, 'config_test.yaml')

            dryrun_return_code = subprocess.call([
                'snakemake',
                '--snakefile', snakefile,
                '--configfile', configfile,
                '-np'
            ])

            self.assertEqual(0, dryrun_return_code)

    def test_chip_pe_dryrun_passes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            os.chdir(temp_dir)

            results_dir = os.path.join(temp_dir, 'results')
            results_temp_dir = os.path.join(results_dir, 'tmp')
            chip_general_config = os.path.join(PIPELINE_BASE_DIR, 'config', 'ChIP_histone_general_pe.yaml')

            ########################################
            # Create the configfile
            config_creator_script = os.path.join(PIPELINE_BASE_DIR, 'scripts', 'config_creator.py')
            general_input = chip_general_config

            test_samplesheet = '''lib,sample,genome,input,homer_fmg_genome,basepath
SRR6214070,GM12878_NKRF,hg19,GM12878_Input,hg19r,{}/data/sra_chip_test_data/SAMN06121956/
SRR5112307,GM12878_Input,hg19,,,{}/data/sra_chip_test_data/SAMN06122176/''' # FIXME: make these samples match new chip test data when we have it
            test_samplesheet = test_samplesheet.format(PIPELINE_BASE_DIR, PIPELINE_BASE_DIR)
            test_sheet_filename = os.path.join(temp_dir, 'samplesheet_test.csv')
            with open(test_sheet_filename, 'w') as test_sheet_filehandle:
                test_sheet_filehandle.write(test_samplesheet)

            outputfile_config = os.path.join(temp_dir, 'config_test.yaml')

            with open(outputfile_config, 'w') as outfile_handle:
                cfg_create_return_code = subprocess.call([
                    config_creator_script,
                    '--general_input', general_input,
                    '--per_lib_input', test_sheet_filename,
                    '--results_dir', results_dir,
                    '--temp_dir', results_temp_dir
                ], stdout=outfile_handle)

                self.assertEqual(0, cfg_create_return_code)

            ########################################
            # Dryrun test
            snakefile = os.path.join(PIPELINE_BASE_DIR, 'ChIPseq_pe.smk')
            configfile = os.path.join(temp_dir, 'config_test.yaml')

            dryrun_return_code = subprocess.call([
                'snakemake',
                '--snakefile', snakefile,
                '--configfile', configfile,
                '-np'
            ])

            self.assertEqual(0, dryrun_return_code)

    def test_atac_endtoend_passes(self):
        # Using specific temporary directory prefix so that worker nodes on GreatLakes cluster
        # can all access the same temp space. If this is not done, TemporaryDirectory creates them
        # in /tmp, which is not shared between nodes
        with tempfile.TemporaryDirectory(prefix = TEST_TMP_DIR) as temp_dir:
            os.chdir(temp_dir)

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
            outputfile_config = os.path.join(temp_dir, 'config_test.yaml')

            test_samplesheet = '''lib,sample,genome,basepath
1,s1_control,mm10,{}/data/atac_test_data/Sample_1
2,s1_treatment,mm10,{}/data/atac_test_data/Sample_2
3,s2_treatment,mm10,{}/data/atac_test_data/Sample_3
4,s3_treatment,mm10,{}/data/atac_test_data/Sample_4'''
            test_samplesheet = test_samplesheet.format(PIPELINE_BASE_DIR, PIPELINE_BASE_DIR, PIPELINE_BASE_DIR, PIPELINE_BASE_DIR)
            test_sheet_filename = os.path.join(temp_dir, 'samplesheet_test.csv')
            with open(test_sheet_filename, 'w') as test_sheet_filehandle:
                test_sheet_filehandle.write(test_samplesheet)

            with open(outputfile_config, 'w') as outfile_handle:
                cfg_create_return_code = subprocess.call([
                    config_creator_script,
                    '--general_input', general_input,
                    '--per_lib_input', test_sheet_filename,
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
