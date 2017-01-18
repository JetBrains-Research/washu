import subprocess
import unittest
import os

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata/bed'
INTERSECT_SH = os.path.dirname(os.path.abspath(__file__)) + '/../bed/intersect.sh'


# 0      100  200    300  400    500  600    700
# A.bed
# |---------| |---------| |---------| |---------|
# B.bed
#            |------------------|              |--|
# C.bed
#  |-| |-|        |-|          |-|                  |-|

class IntersectTest(unittest.TestCase):
    def test_intersect(self):
        ps = subprocess.Popen(['bash', INTERSECT_SH, TEST_DATA + '/A.bed', TEST_DATA + '/B.bed', TEST_DATA + '/C.bed'],
                              stdout=subprocess.PIPE)
        self.assertEqual("""chr1	150	500
""", ps.communicate()[0].decode("utf-8"))

    def test_intersectAB(self):
        ps = subprocess.Popen(['bash', INTERSECT_SH, TEST_DATA + '/A.bed', TEST_DATA + '/B.bed'],
                              stdout=subprocess.PIPE)
        self.assertEqual("""chr1	150	500
chr1	600	750
""", ps.communicate()[0].decode("utf-8"))

    def test_intersectAC(self):
        ps = subprocess.Popen(['bash', INTERSECT_SH, TEST_DATA + '/A.bed', TEST_DATA + '/C.bed'],
                              stdout=subprocess.PIPE)
        self.assertEqual("""chr1	0	100
chr1	200	300
chr1	400	500
""", ps.communicate()[0].decode("utf-8"))

    def test_intersectBC(self):
        ps = subprocess.Popen(['bash', INTERSECT_SH, TEST_DATA + '/B.bed', TEST_DATA + '/C.bed'],
                              stdout=subprocess.PIPE)
        self.assertEqual("""chr1	150	460
""", ps.communicate()[0].decode("utf-8"))



if __name__ == '__main__':
    unittest.main()
