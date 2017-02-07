import subprocess
import unittest
import os

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata/bed'
METAPEAKS_SH = os.path.dirname(os.path.abspath(__file__)) + '/../bed/metapeaks.sh'


# 0      100  200    300  400    500  600    700
# A.bed
# |---------| |---------| |---------| |---------|
# B.bed
#            |------------------|              |--|
# C.bed
#  |-| |-|        |-|          |-|                  |-|

class MetapeaksTest(unittest.TestCase):
    def test_metapeaks(self):
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
