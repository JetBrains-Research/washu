import unittest
from macs2_logs import macs2_logs
import pandas as pd
import os

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata'


class Macs2LogsTest(unittest.TestCase):
    def test_report(self):
        self.assertTrue(pd.read_csv(TEST_DATA + '/macs2/macs2_report_correct.csv', sep=',')
                        .equals(macs2_logs(TEST_DATA + '/macs2')))


if __name__ == '__main__':
    unittest.main()
