import os
import unittest

import pandas as pd

from logs.bowtie_logs import report

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata'


class BowtieLogsTest(unittest.TestCase):
    def test_report(self):
        self.assertEqual(str(pd.read_csv(TEST_DATA + '/bowtie/bowtie_report_correct.csv', sep=',')).strip(),
                         str(report(TEST_DATA + '/bowtie')).strip())


if __name__ == '__main__':
    unittest.main()
