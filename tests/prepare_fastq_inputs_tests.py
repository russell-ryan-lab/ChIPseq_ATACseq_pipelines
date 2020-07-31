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

        expected_result = set([
            os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799398_1.fastq.gz'),
            os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799398_2.fastq.gz'),
            os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799401_1.fastq.gz'),
            os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799401_2.fastq.gz'),
            os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799402_1.fastq.gz'),
            os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799402_2.fastq.gz'),
            os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799447_1.fastq.gz'),
            os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SRR5799447_2.fastq.gz')
        ])

        actual_result = set(prepare_fastq_inputs.get_fastqs(os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data')))

        self.assertEqual(expected_result, actual_result)

    # Test positive case of group_fastqs_by_metadata()
    def test_group_fastqs_by_metadata(self):
        args = Namespace(
            fastq_dir=os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data'),
            metadata_csv=os.path.join(PIPELINE_BASE_DIR, 'data', 'sra_atac_test_data', 'SraRunTable.txt'),
            group_col='BioSample',
            indiv_col='Run',
            strip_regex='_R*[12]_?[0-9]{0,3}\.fastq\.gz$'
        )

        fastqs = prepare_fastq_inputs.get_fastqs(args.fastq_dir)
        metadata = pd.read_csv(args.metadata_csv)

        expected_result = defaultdict(list, {
            'SAMN07312033': [
                'SRR5799398_1.fastq.gz',
                'SRR5799398_2.fastq.gz'
            ],
            'SAMN07312030': [
                'SRR5799401_1.fastq.gz',
                'SRR5799401_2.fastq.gz',
                'SRR5799402_1.fastq.gz',
                'SRR5799402_2.fastq.gz'
            ],
            'SAMN07312104': [
                'SRR5799447_1.fastq.gz',
                'SRR5799447_2.fastq.gz'
            ]
        })

        actual_result = prepare_fastq_inputs.group_fastqs_by_metadata(fastqs, metadata, args.indiv_col, args.group_col, args.strip_regex)

        self.assertEqual(expected_result, actual_result)

    def test_add_readnum_to_filename(self):
        # Test that warning is issued for filenames with _R1_, etc.
        already_R = 'SRR5799402_R2.fastq.gz'
        self.assertWarns(UserWarning, prepare_fastq_inputs.add_readnum_to_filename, already_R)

        # Positive test-case for SRA style filenames
        sra_type_fn = 'SRR5799402_2.fastq.gz'
        expected_result = 'SRR5799402_R2.fastq.gz'
        actual_result = prepare_fastq_inputs.add_readnum_to_filename(sra_type_fn)
        self.assertEqual(expected_result, actual_result)

        # Positive test-case for generic filenames
        generic_fn = 'SRR5799402.fastq.gz'
        expected_result = 'SRR5799402_R1.fastq.gz'
        actual_result = prepare_fastq_inputs.add_readnum_to_filename(generic_fn)
        self.assertEqual(expected_result, actual_result)

    def test_check_problematic_readnames(self):
        problem_file = os.path.join(PIPELINE_BASE_DIR, 'tests', 'test_files', 'SRR2932619_1.fastq.gz'),
        self.assertWarns(UserWarning, prepare_fastq_inputs.check_problematic_readnames(problem_file))
