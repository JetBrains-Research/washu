import os
import unittest

import pandas as pd

from logs.bowtie_logs import process_bowtie_logs

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata/bowtie'


class BowtieLogsTest(unittest.TestCase):
    def test_process_bowtie_logs(self):
        file = TEST_DATA + '/bowtie_report.csv'
        if os.path.exists(file):
            os.remove(file)
        try:
            process_bowtie_logs(TEST_DATA)
            self.assertTrue(os.path.exists(file))
            self.assertEqual(str(pd.read_csv(TEST_DATA + '/bowtie_report_correct.csv', sep=',')).strip(),
                             str(pd.read_csv(file, sep=',')).strip())
        finally:
            os.remove(file)


if __name__ == '__main__':
    unittest.main()
