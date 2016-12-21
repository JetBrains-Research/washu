import os
import unittest

import pandas as pd

from logs.bowtie2_logs import bowtie2_logs

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata'


class Bowtie2LogsTest(unittest.TestCase):
    def test_report(self):
        self.assertEqual(str(pd.read_csv(TEST_DATA + '/bowtie2/bowtie2_report_correct.csv', sep=',')),
                         str(bowtie2_logs(TEST_DATA + '/bowtie2')))


if __name__ == '__main__':
    unittest.main()
