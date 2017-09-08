import os
from pathlib import Path
from bed.bedtrace import Bed, union, intersect, minus, compare

import pytest
from test.fixtures import test_data, bedtrace_cleanup


# 0      100  200    300  400    500  600    700
# A.bed
# |---------| |---------| |---------| |---------|
# B.bed
#            |------------------|              |--|
# C.bed
#  |-| |-|        |-|          |-|                  |-|
def test_str_simple(test_data):
    assert str(Bed(test_data("bed/A.bed"))) == "A.bed"


def test_str_compound(test_data):
    c = minus(union(Bed(test_data("bed/B.bed")),
                    Bed(test_data("bed/A.bed"))),
              Bed(test_data("bed/C.bed")))
    assert """minus
\tunion
\t\tB.bed
\t\tA.bed
\tC.bed""" == str(c)


def test_union(test_data):
    u = union(Bed(test_data("bed/A.bed")),
              Bed(test_data("bed/B.bed")))
    assert u.path is None

    u.compute()
    assert u.path is not None

    actual = Path(u.path).read_text().replace(test_data("bed/"), '')
    assert actual == """chr1	0	100	1_A.bed
chr1	150	500	1_A.bed|2_B.bed
chr1	600	750	1_A.bed|2_B.bed
"""


def test_intersect(test_data):
    i = intersect(Bed(test_data("bed/A.bed")), Bed(test_data("bed/B.bed")))
    assert i.path is None

    i.compute()
    assert i.path is not None
    assert Path(i.path).read_text() == ("chr1	150	500\n"
                                        "chr1	600	750\n")


def test_minus(test_data):
    m = minus(Bed(test_data("bed/A.bed")), Bed(test_data("bed/B.bed")))
    assert m.path is None
    m.compute()
    assert m.path is not None
    assert Path(m.path).read_text() == "chr1	0	100\n"


@pytest.mark.parametrize("files,expected_count", [
    (["A.bed"], 4),
    (["B.bed"], 2),
    (["C.bed"], 5),
    (["A.bed", "B.bed"], 3),
])
def test_count(test_data, files, expected_count):
    if len(files) == 1:
        bed = Bed(test_data("bed/" + files[0]))
    else:
        bed = union(*[Bed(test_data("bed/" + f)) for f in files])
    assert bed.count() == expected_count


def test_compare(test_data):
    c = compare(Bed(test_data("bed/A.bed")), Bed(test_data("bed/B.bed")))
    assert c.path is None

    c.compute()
    assert c.path is not None
    assert Path(c.cond1).read_text() == "chr1	0	100\n"
    assert Path(c.cond2).read_text() == ""
    assert Path(c.common).read_text() == ("chr1	150	500\n"
                                          "chr1	600	750\n")


def test_save(test_data):
    assert union(Bed(test_data("bed/A.bed")),
                 Bed(test_data("bed/B.bed"))).count() == 3


def test_process_pvalue(test_data):
    u = union(Bed(test_data("bed/A.narrowPeak")),
              Bed(test_data("bed/B.narrowPeak")))

    f = u.process_pvalue()
    assert Path(f).read_text() == """chr2	9745187	9746077	26
chr2	9746391	9746765	20
chr1	4857963	4858364	6
chr1	4807879	4808181	4
"""


def test_process_pvalue_single(test_data):
    u = minus(Bed(test_data("bed/A.narrowPeak")),
              intersect(Bed(test_data("bed/A.narrowPeak")),
                        Bed(test_data("bed/B.narrowPeak"))))
    f = u.process_pvalue()
    assert Path(f).read_text() == """chr2	9745187	9746077	26
chr2	9746391	9746765	20
"""


def test_process_diffbind(test_data):
    u = minus(Bed(test_data("bed/A_diffbind.bed")),
              Bed(test_data("bed/A.narrowPeak")))
    f = u.process_pvalue()
    assert Path(f).read_text() == """chr1	28417357	28419402	1.53e-06
chr1	8212454	8213531	7.98e-10
"""


def test_process_bdgdiff(test_data):
    u = Bed(test_data("bed/A_B_bdgdiff_cond1.bed"))
    f = u.process_pvalue()
    assert Path(f).read_text() == """chr1	92780966	92781404	5.47439
chr1	10571239	10571874	4.94096
"""


def test_save3(test_data):
    u = Bed(test_data("bed/A_B_bdgdiff_cond1.bed"))
    # TODO test save3
    f = u.process_pvalue()
    assert Path(f).read_text() == """chr1	92780966	92781404	5.47439
chr1	10571239	10571874	4.94096
"""
