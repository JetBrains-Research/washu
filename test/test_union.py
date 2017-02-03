import subprocess
import unittest
import os

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata/bed'
UNION_SH = os.path.dirname(os.path.abspath(__file__)) + '/../bed/union.sh'


# 0      100  200    300  400    500  600    700
# A.bed
# |---------| |---------| |---------| |---------|
# B.bed
#            |------------------|              |--|
# C.bed
#  |-| |-|        |-|          |-|                  |-|

class UnionTest(unittest.TestCase):
    def test_union(self):
        ps = subprocess.Popen(['bash', UNION_SH, TEST_DATA + '/A.bed', TEST_DATA + '/B.bed', TEST_DATA + '/C.bed'],
                              stdout=subprocess.PIPE)
        self.assertEqual("""chr1	0	100	1|3|3
chr1	150	500	2|1|3|1|3
chr1	600	750	1|2
chr1	800	850	3
""", ps.communicate()[0].decode("utf-8"))

if __name__ == '__main__':
    unittest.main()
