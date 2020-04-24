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

    def test_read_input(self):
        invalid_yaml = os.path.join(PIPELINE_BASE_DIR, 'tests', 'test_files', 'invalid.yaml')
        invalid_json = os.path.join(PIPELINE_BASE_DIR, 'tests', 'test_files', 'invalid.json')

        self.assertRaisesRegex(SystemExit, 'Error loading YAML', config_creator.read_input, invalid_yaml)
        self.assertRaisesRegex(SystemExit, 'Error loading JSON', config_creator.read_input, invalid_json)
        self.assertRaises(FileNotFoundError, config_creator.read_input, 'foo.yaml')
        self.assertRaises(FileNotFoundError, config_creator.read_input, 'foo.json')
