import tempfile
import unittest

import pandas as pd
import scripts.config_creator as config_creator


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

    def test_read_input(self):
        invalid_yaml = '@: foo'
        invalid_json = 'foo: bar'

        loaded_input = config_creator.read_input()
