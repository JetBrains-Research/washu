import subprocess
import unittest
import os

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata/bed'
MINUS_SH = os.path.dirname(os.path.abspath(__file__)) + '/../bed/minus.sh'


# 0      100  200    300  400    500  600    700
# A.bed
# |---------| |---------| |---------| |---------|
# B.bed
#            |------------------|              |--|
# C.bed
#  |-| |-|        |-|          |-|                  |-|

class MinusTest(unittest.TestCase):
    def test_minusAB(self):
        ps = subprocess.Popen(['bash', MINUS_SH, TEST_DATA + '/A.bed', TEST_DATA + '/B.bed'],
                              stdout=subprocess.PIPE)
        self.assertEqual("""chr1	0	100
""", ps.communicate()[0].decode("utf-8"))

    def test_minusAC(self):
        ps = subprocess.Popen(['bash', MINUS_SH, TEST_DATA + '/A.bed', TEST_DATA + '/C.bed'],
                              stdout=subprocess.PIPE)
        self.assertEqual("""chr1	600	700
""", ps.communicate()[0].decode("utf-8"))

    def test_minusBC(self):
        ps = subprocess.Popen(['bash', MINUS_SH, TEST_DATA + '/B.bed', TEST_DATA + '/C.bed'],
                              stdout=subprocess.PIPE)
        self.assertEqual("""chr1	650	750
""", ps.communicate()[0].decode("utf-8"))


if __name__ == '__main__':
    unittest.main()
