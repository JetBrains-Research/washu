import os
import unittest
from pathlib import Path

from bed.metabed import Bed, union, intersect, minus, compare, cleanup

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata/bed'
UNION_SH = os.path.dirname(os.path.abspath(__file__)) + '/../bed/union.sh'


# 0      100  200    300  400    500  600    700
# A.bed
# |---------| |---------| |---------| |---------|
# B.bed
#            |------------------|              |--|
# C.bed
#  |-| |-|        |-|          |-|                  |-|

class MetaBedTest(unittest.TestCase):
    def test_str(self):
        self.assertEqual("A.bed", str(Bed(TEST_DATA + '/A.bed')))
        self.assertEqual("union(A.bed, B.bed)", str(union(Bed(TEST_DATA + '/A.bed'), Bed(TEST_DATA + '/B.bed'))))
        self.assertEqual("union(A.bed, B.bed, C.bed)",
                         str(union(Bed(TEST_DATA + '/B.bed'), Bed(TEST_DATA + '/A.bed'), Bed(TEST_DATA + '/C.bed'))))

    def test_union(self):
        u = union(Bed(TEST_DATA + '/A.bed'), Bed(TEST_DATA + '/B.bed'))
        self.assertIsNone(u.path)
        u.compute()
        self.assertIsNotNone(u.path)
        self.assertEqual("""chr1	0	100	1
chr1	150	500	2|1|1
chr1	600	750	1|2
""", Path(u.path).read_text())

    def test_intersect(self):
        i = intersect(Bed(TEST_DATA + '/A.bed'), Bed(TEST_DATA + '/B.bed'))
        self.assertIsNone(i.path)
        i.compute()
        self.assertIsNotNone(i.path)
        self.assertEqual("""chr1	150	500
chr1	600	750
""", Path(i.path).read_text())

    def test_minus(self):
        m = minus(Bed(TEST_DATA + '/A.bed'), Bed(TEST_DATA + '/B.bed'))
        self.assertIsNone(m.path)
        m.compute()
        self.assertIsNotNone(m.path)
        self.assertEqual("""chr1	0	100
""", Path(m.path).read_text())

    def test_count(self):
        self.assertEqual(4, Bed(TEST_DATA + '/A.bed').count())
        self.assertEqual(2, Bed(TEST_DATA + '/B.bed').count())
        self.assertEqual(5, Bed(TEST_DATA + '/C.bed').count())
        self.assertEqual(3, union(Bed(TEST_DATA + '/A.bed'), Bed(TEST_DATA + '/B.bed')).count())

    def test_compare(self):
        c = compare(Bed(TEST_DATA + '/A.bed'), Bed(TEST_DATA + '/B.bed'))
        self.assertIsNone(c.path)
        c.compute()
        self.assertIsNotNone(c.path)
        self.assertEqual("""chr1	0	100
""", Path(c.cond1).read_text())
        self.assertEqual('', Path(c.cond2).read_text())
        self.assertEqual("""chr1	150	500
chr1	600	750
""", Path(c.common).read_text())

    def test_save(self):
        self.assertEqual(3, union(Bed(TEST_DATA + '/A.bed'), Bed(TEST_DATA + '/B.bed')).count())

    @classmethod
    def tearDownClass(cls):
        cleanup()


if __name__ == '__main__':
    unittest.main()
