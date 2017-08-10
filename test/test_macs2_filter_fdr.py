import os
import shutil
import unittest

from bed.bedtrace import run

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata/macs2_filter_fdr'
MACS2_FILTER_FDR_SH = os.path.dirname(os.path.abspath(__file__)) + '/../bed/macs2_filter_fdr.sh'


class TestMacs2FilterFdr(unittest.TestCase):
    def testMacs2FilterFdr(self):
        try:
            run([['bash', MACS2_FILTER_FDR_SH, TEST_DATA, TEST_DATA + "/test", "0.1", "0.01"]])
            for dirpath, dirs, files in os.walk(TEST_DATA + "/test"):
                peaks = [f for f in files if "0.01" in f]
                self.assertTrue(len(peaks) > 0)
                self.assertTrue(peaks[0].endswith('.broadPeak'))
                with open(TEST_DATA + "/test/" + peaks[0]) as peak:
                    line = [l for l in peak][1]
                    self.assertEqual("chr1\t241186\t241332\tfoo_0.01_peak_2\t20\t.\t2.82218\t3.36676\t2.05718\t\n", line)
            run([['bash', MACS2_FILTER_FDR_SH, TEST_DATA, TEST_DATA + "/test", "0.1", "0.05"]])
            for dirpath, dirs, files in os.walk(TEST_DATA + "/test"):
                peaks = [f for f in files if "0.05" in f]
                self.assertTrue(len(peaks) > 0)
                self.assertTrue(peaks[0].endswith('.broadPeak'))
                with open(TEST_DATA + "/test/" + peaks[0]) as peak:
                    line = [l for l in peak][1]
                    self.assertEqual("chr1\t241186\t241332\tfoo_0.05_peak_2\t20\t.\t2.82218\t3.36676\t2.05718\t\n", line)
        finally:
            pass
            shutil.rmtree(TEST_DATA + "/test")


if __name__ == '__main__':
    unittest.main()
