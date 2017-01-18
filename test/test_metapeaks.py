import subprocess
import unittest
import os

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata/bed'
METAPEAKS_SH = os.path.dirname(os.path.abspath(__file__)) + '/../bed/metapeaks.sh'


class MetapeaksTest(unittest.TestCase):
    def test_metapeaks2(self):
        ps = subprocess.Popen(['bash', METAPEAKS_SH, TEST_DATA + '/A.bed', TEST_DATA + '/B.bed'],
                              stdout=subprocess.PIPE)
        self.assertEqual("""PEAKS:        4        2
1 0	1
1 1	2
1, 2
""", ps.communicate()[0].decode("utf-8"))

    def test_metapeaks3(self):
        ps = subprocess.Popen(['bash', METAPEAKS_SH, TEST_DATA + '/A.bed', TEST_DATA + '/B.bed', TEST_DATA + '/C.bed'],
                              stdout=subprocess.PIPE)
        self.assertEqual("""PEAKS:        4        2        5
0 0 1	1
1 0 1	1
1 1 0	1
1 1 1	1
1, 1, 1, 1
""", ps.communicate()[0].decode("utf-8"))

if __name__ == '__main__':
    unittest.main()
