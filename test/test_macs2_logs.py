import os
import unittest

import pandas as pd

from logs.macs2_logs import report

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata'


class Macs2LogsTest(unittest.TestCase):
    def test_report(self):
        self.assertEqual(str(pd.read_csv(TEST_DATA + '/macs2/macs2_report_correct.csv', sep=',')).strip(),
                         str(report(TEST_DATA + '/macs2')).strip())


if __name__ == '__main__':
    unittest.main()
