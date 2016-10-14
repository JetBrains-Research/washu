import unittest
from bowtie_logs import bowtie_logs
import pandas as pd
import os

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata'


class BowtieLogsTest(unittest.TestCase):
    def test_report(self):
        self.assertEqual(str(pd.read_csv(TEST_DATA + '/bowtie/bowtie_report_correct.csv', sep=',')),
                         str(bowtie_logs(TEST_DATA + '/bowtie')))


if __name__ == '__main__':
    unittest.main()
