from pathlib import Path

import pytest

from bed.bedtrace import Bed, union, intersect, minus, compare, jaccard
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


def test_jaccard(test_data):
    u = jaccard(test_data("bed/A.bed"), test_data("bed/B.bed"))
    assert u == 1.0 / 3.0


def test_jaccard_not_sorted(test_data):
    u = jaccard(test_data("bed/A.unsorted.bed"), test_data("bed/B.unsorted.bed"))
    assert u == 1.0 / 3.0


def test_jaccard_self_intersections(test_data):
    u = jaccard(test_data("bed/E.bed"), test_data("bed/F.bed"))
    assert u == 35.0 / 72.0


def test_save(test_data):
    assert union(Bed(test_data("bed/A.bed")),
                 Bed(test_data("bed/B.bed"))).count() == 3


def test_cat(test_data):
    u = Bed(test_data("bed/A.bed"))

    assert u.cat() == """chr1	0	100
chr1	200	300
chr1	400	500
chr1	600	700
"""
