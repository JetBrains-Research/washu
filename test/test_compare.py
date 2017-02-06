import subprocess
import tempfile
import unittest
import os

from pathlib import Path

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata/bed'
COMPARE_SH = os.path.dirname(os.path.abspath(__file__)) + '/../bed/compare.sh'


# 0      100  200    300  400    500  600    700
# A.bed
# |---------| |---------| |---------| |---------|
# B.bed
#            |------------------|              |--|
# C.bed
#  |-| |-|        |-|          |-|                  |-|

class CompareTest(unittest.TestCase):
    def test_compareAB(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', prefix='metabeds', delete=False) as tmpfile:
            prefix = tmpfile.name.replace('.txt', '')
            subprocess.Popen(['bash', COMPARE_SH, TEST_DATA + '/A.bed', TEST_DATA + '/B.bed', prefix],
                             stdout=subprocess.PIPE).communicate()
            cond1 = prefix + "_cond1.bed"
            cond2 = prefix + "_cond2.bed"
            common = prefix + "_common.bed"
        self.assertEqual("""chr1	0	100
""", Path(cond1).read_text())
        self.assertEqual('', Path(cond2).read_text())
        self.assertEqual("""chr1	150	500
chr1	600	750
""", Path(common).read_text())


if __name__ == '__main__':
    unittest.main()
