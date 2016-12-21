import os
import unittest

import pandas as pd

from logs.macs2_logs import macs2_logs

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata'


class Macs2LogsTest(unittest.TestCase):
    def test_report(self):
        self.assertEqual(str(pd.read_csv(TEST_DATA + '/macs2/macs2_report_correct.csv', sep=',')),
                         str(macs2_logs(TEST_DATA + '/macs2')))


if __name__ == '__main__':
    unittest.main()
