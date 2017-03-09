import os
import unittest
from pathlib import Path

from bed.bedtrace import Bed, union, intersect, minus, compare, cleanup

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata/bed'
UNION_SH = os.path.dirname(os.path.abspath(__file__)) + '/../bed/union.sh'


# 0      100  200    300  400    500  600    700
# A.bed
# |---------| |---------| |---------| |---------|
# B.bed
#            |------------------|              |--|
# C.bed
#  |-| |-|        |-|          |-|                  |-|

class BedTrace(unittest.TestCase):
    def test_str(self):
        self.assertEqual("A.bed", str(Bed(TEST_DATA + '/A.bed')))
        self.assertEqual("""minus
	union
		B.bed
		A.bed
	C.bed""", str(minus(
            union(Bed(TEST_DATA + '/B.bed'), Bed(TEST_DATA + '/A.bed')),
            Bed(TEST_DATA + '/C.bed'))))

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

    def test_process_pvalue(self):
        u = union(Bed(TEST_DATA + '/A.narrowPeak'), Bed(TEST_DATA + '/B.narrowPeak'))
        f = u.process_pvalue()
        self.assertEqual("""chr2	9745187	9746077	26
chr2	9746391	9746765	20
chr1	4857963	4858364	6
chr1	4807879	4808181	4
""", Path(f).read_text())

    def test_process_pvalue_single(self):
        u = minus(Bed(TEST_DATA + '/A.narrowPeak'),
                  intersect(Bed(TEST_DATA + '/A.narrowPeak'), Bed(TEST_DATA + '/B.narrowPeak')))
        f = u.process_pvalue()
        self.assertEqual("""chr2	9745187	9746077	26
chr2	9746391	9746765	20
""", Path(f).read_text())

    def test_process_diffbind(self):
        u = minus(Bed(TEST_DATA + '/A_diffbind.bed'), Bed(TEST_DATA + '/A.narrowPeak'))
        f = u.process_pvalue()
        self.assertEqual("""chr1	28417357	28419402	1.53e-06
chr1	8212454	8213531	7.98e-10
""", Path(f).read_text())

    def test_process_bdgdiff(self):
        u = Bed(TEST_DATA + '/A_B_bdgdiff_cond1.bed')
        f = u.process_pvalue()
        self.assertEqual("""chr1	92780966	92781404	5.47439
chr1	10571239	10571874	4.94096
""", Path(f).read_text())

    def test_save3(self):
        u = Bed(TEST_DATA + '/A_B_bdgdiff_cond1.bed')
        # TODO test save3
        f = u.process_pvalue()
        self.assertEqual("""chr1	92780966	92781404	5.47439
chr1	10571239	10571874	4.94096
""", Path(f).read_text())


    @classmethod
    def tearDownClass(cls):
        pass
        # cleanup()


if __name__ == '__main__':
    unittest.main()
