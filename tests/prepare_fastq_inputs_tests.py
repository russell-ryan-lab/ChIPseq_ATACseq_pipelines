from collections import defaultdict
import json
import os
import tempfile
import unittest
import yaml

import pandas as pd
from argparse import Namespace

import scripts.prepare_fastq_inputs as prepare_fastq_inputs

TEST_DIR = os.path.realpath(os.path.dirname(__file__))
PIPELINE_BASE_DIR = os.path.abspath(os.path.join(TEST_DIR, '..'))

class PrepareFastqInputsTests(unittest.TestCase):
    def setUp(self):
        self.original_wd = os.getcwd()

    def tearDown(self):
        os.chdir(self.original_wd)

    def test_validate_input_arguments(self):
        args = Namespace(
            fastq_dir=os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data'),
            metadata_csv=os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SraRunTable.txt'),
            group_col='BioSample',
            indiv_col='Run'
        )

        metadata = pd.read_csv(args.metadata_csv)

        args_no_fastqs = Namespace(
            fastq_dir=os.path.join(PIPELINE_BASE_DIR, 'blue'),
            metadata_csv=os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SraRunTable.txt'),
            group_col='BioSample',
            indiv_col='Run'
        )

        metadata_sans_Run = metadata.drop(['Run'], axis = 1)
        metadata_sans_BioSample = metadata.drop(['BioSample'], axis = 1)

        # Test for missing metadata columns (validate_input_arguments)
        self.assertRaises(RuntimeError, prepare_fastq_inputs.validate_input_arguments, metadata_sans_Run, args)
        self.assertRaises(RuntimeError, prepare_fastq_inputs.validate_input_arguments, metadata_sans_BioSample, args)

        # Test for valid fastq directory (validate_input_arguments)
        self.assertRaises(RuntimeError, prepare_fastq_inputs.validate_input_arguments, metadata, args_no_fastqs)

    # Test for error thrown with no fastqs in dir
    def test_get_fastqs(self):
        self.assertRaises(RuntimeError, prepare_fastq_inputs.get_fastqs, os.path.join(PIPELINE_BASE_DIR, 'config'))

        expected_result = set([os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799398_1.fastq'), os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799398_2.fastq'), os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799401_1.fastq'), os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799401_2.fastq'), os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799402_1.fastq'), os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799402_2.fastq'), os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799447_1.fastq'), os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799447_2.fastq')])

        actual_result = set(prepare_fastq_inputs.get_fastqs(os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data')))

        self.assertEqual(expected_result, actual_result)

    # Test positive case of group_fastqs_by_metadata()
    def test_group_fastqs_by_metadata(self):
        args = Namespace(
            fastq_dir=os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data'),
            metadata_csv=os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SraRunTable.txt'),
            group_col='BioSample',
            indiv_col='Run',
            strip_regex='_R*[12](\.fastq|\.fq)(\.gz)*$'
        )

        fastqs = prepare_fastq_inputs.get_fastqs(args.fastq_dir)
        metadata = pd.read_csv(args.metadata_csv)

        expected_result = defaultdict(list, {'SAMN07312033': ['SRR5799398_2.fastq', 'SRR5799398_1.fastq'], 'SAMN07312030': ['SRR5799401_1.fastq', 'SRR5799401_2.fastq', 'SRR5799402_2.fastq', 'SRR5799402_1.fastq'], 'SAMN07312104': ['SRR5799447_1.fastq', 'SRR5799447_2.fastq']})

        actual_result = prepare_fastq_inputs.group_fastqs_by_metadata(fastqs, metadata, args.indiv_col, args.group_col, args.strip_regex)

        self.assertEqual(expected_result, actual_result)
