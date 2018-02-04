import pytest
from pathlib import Path
from downstream.signals import signals_util as su


@pytest.mark.parametrize("file,norm", [
    ("/foo/H3K27ac/hg19_100000/hg19_100000_diffbind_tmm_reads_effective_cpm.tsv",
     "diffbind_tmm_reads_effective_cpm"),
    ("/foo/H3K27ac/hg19_100000/hg19_100000_rawz.tsv", "rawz"),
    ("/foo/k27me3@dmrs/k27me3@dmrs_dedup_FALSE_f_125_DBA_SCORE_READS_counts.csv",
     "dedup_FALSE_f_125_DBA_SCORE_READS_counts"),
    ("/foo/hg19_10000/boo.doo.csv", "boo.doo")
])
def test_extract_normalization(file, norm):
    assert su.extract_normalization(Path(file)) == norm


@pytest.mark.parametrize("file,dtype", [
    ("/foo/H3K27ac/hg19_100000/foo.tsv", "H3K27ac"),
    ("/foo/meth/hg19_100000/foo.tsv", "meth"),
    ("/foo/boo/h3k27me3/hg19_100000/foo.tsv", "h3k27me3"),
    ("/foo/boo/k27me3/hg19_100000/foo.tsv", "N/A"),
    ("/foo/boo/hg19_100000/foo.tsv", "N/A"),
])
def test_extract_datatype(file, dtype):
    assert su.extract_datatype(Path(file)) == dtype
