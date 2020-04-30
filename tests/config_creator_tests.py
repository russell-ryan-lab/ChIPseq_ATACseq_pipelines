import json
import os
import tempfile
import unittest
import yaml

import pandas as pd
try:
    #pylint: disable=locally-disabled,unused-import
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import scripts.config_creator as config_creator

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
PIPELINE_BASE_DIR = os.path.abspath(os.path.join(TEST_DIR, '..'))


class ConfigCreatorTests(unittest.TestCase):
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

    def test_validate_sample_table(self):
        dat_missingcol = {'sample': pd.Series(['s1_foo', 's1_2*(@%', 's2_bar', 's3_biz']),
        'basepath' : pd.Series(['foo', 'bar', 'baz', 'bin'])
        }
        sample_table_missingcol = pd.DataFrame(dat_missingcol)
        self.assertRaises(RuntimeError, config_creator.validate_sample_table, sample_table_missingcol)

        dat_badname = {'sample': pd.Series(['s1_foo', 's1_2*(@%', 's2_bar', 's3_biz']),
        'lib' : pd.Series(['1234', '2345', '3456', '4567']),
        'basepath' : pd.Series(['foo', 'bar', 'baz', 'bin'])
        }
        sample_table_badname = pd.DataFrame(dat_badname)
        self.assertRaises(RuntimeError, config_creator.validate_sample_table, sample_table_badname)

    def test_read_input(self):
        invalid_yaml = os.path.join(PIPELINE_BASE_DIR, 'tests', 'test_files', 'invalid.yaml')
        invalid_json = os.path.join(PIPELINE_BASE_DIR, 'tests', 'test_files', 'invalid.json')

        valid_yaml = os.path.join(PIPELINE_BASE_DIR, 'config', 'ATAC_general.yaml')
        valid_json = os.path.join(PIPELINE_BASE_DIR, 'config', 'ATAC_general.json')

        # Test error-handling
        self.assertRaisesRegex(SystemExit, 'Error loading YAML', config_creator.read_input, invalid_yaml)
        self.assertRaisesRegex(SystemExit, 'Error loading JSON', config_creator.read_input, invalid_json)
        self.assertRaises(FileNotFoundError, config_creator.read_input, 'foo.yaml')
        self.assertRaises(FileNotFoundError, config_creator.read_input, 'foo.json')

    def test_basepath_to_filepathsdict(self):
        default_glob = '*.fastq.gz'
        default_capture_regex = '.*_R(\d+).*\.fastq\.gz'
        sample_1_fastqs_dir = os.path.join(PIPELINE_BASE_DIR, 'data', 'atac_test_data', 'Sample_1')
        sample_2_fastqs_dir = os.path.join(PIPELINE_BASE_DIR, 'data', 'atac_test_data', 'Sample_2')

        no_fastqs_dir = os.path.join(PIPELINE_BASE_DIR, 'scripts')
        wrong_glob = 'blue'
        wrong_capture_regex = 'blue'

        self.assertRaises(RuntimeError, config_creator.basepath_to_filepathsdict, no_fastqs_dir, default_glob, default_capture_regex)
        self.assertRaises(RuntimeError, config_creator.basepath_to_filepathsdict, sample_1_fastqs_dir, wrong_glob, default_capture_regex)
        self.assertRaises(RuntimeError, config_creator.basepath_to_filepathsdict, sample_1_fastqs_dir, default_glob, wrong_capture_regex)

        expected_result_1 = {'1': ['/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/ChIPseq_ATACseq_pipelines/data/atac_test_data/Sample_1/1_L001_R1.fastq.gz'], '2': ['/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/ChIPseq_ATACseq_pipelines/data/atac_test_data/Sample_1/1_L001_R2.fastq.gz']}

        actual_result_1 = config_creator.basepath_to_filepathsdict(sample_1_fastqs_dir, default_glob, default_capture_regex)

        self.assertEqual(expected_result_1, actual_result_1)

        expected_result_2 = {'1': ['/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/ChIPseq_ATACseq_pipelines/data/atac_test_data/Sample_2/2_L001_R1.fastq.gz', '/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/ChIPseq_ATACseq_pipelines/data/atac_test_data/Sample_2/2_L002_R1.fastq.gz'], '2': ['/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/ChIPseq_ATACseq_pipelines/data/atac_test_data/Sample_2/2_L001_R2.fastq.gz', '/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/ChIPseq_ATACseq_pipelines/data/atac_test_data/Sample_2/2_L002_R2.fastq.gz']}

        actual_result_2 = config_creator.basepath_to_filepathsdict(sample_2_fastqs_dir, default_glob, default_capture_regex)

        self.assertEqual(expected_result_2, actual_result_2)
