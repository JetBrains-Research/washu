from pathlib import Path
import pytest
from test.fixtures import test_data, tmp_dir
import reports.loci_of_interest as loi


def test_collect_chromhmm(tmp_dir):
    loci_root = Path(tmp_dir)
    generate_test_data_chromhmm(loci_root)
    (loci_root / "foo.bed").touch()
    for folder in ["enhancers", "tfs", "regulatory", "repeats", 'chromhmm']:
        (loci_root / folder).mkdir(parents=True, exist_ok=True)

    res = loi.collect_loci(loci_root)["chromhmm"]
    assert 5 == len(res)
    expected_fnames = [
        'cd14_chromhmm.hg19.1_TssA.bed', 'cd14_chromhmm.hg19.2_TssFlnk.bed',
        'cd14_chromhmm.hg19.9_EnhA1.bed', 'cd14_chromhmm.hg19.10_EnhA2.bed',
        'cd14_chromhmm.hg19.12_ZNF_Rpts.bed'
    ]
    assert expected_fnames == [f.name for f in res]


@pytest.mark.parametrize("fname,descr", [
    ("cd14_chromhmm.hg19.1_TssA.bed", "1_TssA (Active TSS)"),
    ("cd14_chromhmm.hg19.12_ZNF_Rpts.bed", "12_ZNF_Rpts (ZNF genes & repeats)"),
    ("boo.bed", "boo.bed"),
])
def test_chromhmm_state_descr(fname, descr):
    assert descr == loi.chromhmm_state_descr(fname)


def test_collect_loci(tmp_dir):
    loci_root = Path(tmp_dir)
    generate_test_data_chromhmm(loci_root)
    (loci_root / "foo.bed").touch()
    (loci_root / "doo.bed").touch()
    (loci_root / "aaa").touch()
    (loci_root / "aaa.csv").touch()
    (loci_root / "chipseq_diff_loci/H3K27ac").mkdir(parents=True)
    (loci_root / "chipseq_diff_loci/H3K27ac/diff1.bed").touch()
    (loci_root / "chipseq_diff_loci/H3K27ac/diff2.bed").touch()
    for i, name in enumerate(["enhancers", "tfs", "regulatory", "repeats",
                              "golden_median_consensus", "zinbra_median_consensus",
                              "golden_consensus", "zinbra_consensus", 'else']):
        folder = loci_root / name
        folder.mkdir()
        (folder / "boo.bed").touch()
        (folder / "{}.bed".format(i + 1)).touch()

    table = loi.collect_loci(loci_root)
    assert 13 == len(table)
    assert ['None', 'chipseq_diff_loci', 'chromhmm', 'default', 'else', 'enhancers',
            'golden_consensus', 'golden_median_consensus',
            'regulatory', 'repeats', 'tfs',
            'zinbra_consensus', 'zinbra_median_consensus'
            ] == sorted(str(k) for k in table)
    assert 27 == len(table[None])
    assert 2 == len(table['chipseq_diff_loci'])
    assert 2 == len(table['zinbra_consensus'])
    assert 2 == len(table['zinbra_median_consensus'])
    assert ['1.bed', '2.bed', '3.bed', '4.bed', '5.bed', '6.bed', '7.bed', '8.bed', '9.bed',
            'boo.bed', 'boo.bed', 'boo.bed', 'boo.bed', 'boo.bed', 'boo.bed', 'boo.bed',
            'boo.bed', 'boo.bed',
            'cd14_chromhmm.hg19.10_EnhA2.bed', 'cd14_chromhmm.hg19.12_ZNF_Rpts.bed',
            'cd14_chromhmm.hg19.1_TssA.bed', 'cd14_chromhmm.hg19.2_TssFlnk.bed',
            'cd14_chromhmm.hg19.9_EnhA1.bed',
            'diff1.bed', 'diff2.bed', 'doo.bed', 'foo.bed'] == (
        [t.name for t in table[None]]
    )
    assert 15 == len(table['default'])
    assert ['1.bed', '2.bed', '3.bed', '4.bed', 'boo.bed', 'boo.bed', 'boo.bed', 'boo.bed',
            'cd14_chromhmm.hg19.10_EnhA2.bed', 'cd14_chromhmm.hg19.12_ZNF_Rpts.bed',
            'cd14_chromhmm.hg19.1_TssA.bed', 'cd14_chromhmm.hg19.2_TssFlnk.bed',
            'cd14_chromhmm.hg19.9_EnhA1.bed',
            'doo.bed', 'foo.bed'] == [t.name for t in table['default']]


@pytest.mark.parametrize("root,relative_path,selected,mod", [
    # Zinbra
    ("peaks/H3K27ac", "YD_YD4_H3K27ac_hg19_1.0E-6_peaks.bed", True, "H3K27ac"),
    ("peaks/H3K27ac", "OD_OD7_H3K27ac_hg19_1.0E-6_peaks.bed", True, "H3K27ac"),
    ("peaks/H3K27ac", "OD_OD7_H3K27ac_hg19_1.0E-6_peaks.bed_rip.csv", False, "H3K27ac"),
    ("peaks/H3K27ac", "H3K27ac_zinbra_ODS_weak_consensus.bed", False, "H3K27ac"),
    ("peaks/H3K27ac", "H3K27ac_zinbra_weak_consensus.bed", False, "H3K27ac"),
    ("peaks/H3K27ac", "H3K27ac_zinbra_weak_consensus.bed", False, "H3K27ac"),
])
def test_collect_zinbra_peaks(tmp_dir, root, relative_path, selected, mod):
    data_root = Path(tmp_dir)
    file = data_root / root / relative_path
    file.parent.mkdir(parents=True)
    file.touch()

    peaks = loi._collect_zinbra_peaks(data_root)
    assert selected == (file in peaks[mod])


@pytest.mark.parametrize("root,relative_path,selected,mod,exclude_outliers", [
    # Macs2 broad: new layout
    ("H3K27ac", "YD9_k27ac_hg19_broad_peaks.broadPeak", True, "H3K27ac", True),
    ("H3K27ac", "OD9_k27ac_hg19_broad_peaks.broadPeak", True, "H3K27ac", True),
    ("H3K27ac", "YD8_k27ac_hg19_broad_peaks.broadPeak_rip.csv", False, "H3K27ac", True),

    # Macs2 broad
    ("H3K27ac/bed", "YD9_k27ac_hg19_broad_peaks.broadPeak", True, "H3K27ac", True),
    ("H3K27ac/bed", "OD9_k27ac_hg19_broad_peaks.broadPeak", True, "H3K27ac", True),
    ("H3K27ac/bed", "YD8_k27ac_hg19_broad_peaks.broadPeak_rip.csv", False, "H3K27ac", True),

    # Macs2 narrow
    ("H3K4me3/bed", "YD4_k4me3_hg19_fdr_peaks.narrowPeak", True, "H3K4me3", True),
    ("H3K4me3/bed", "OD4_k4me3_hg19_fdr_peaks.narrowPeak", True, "H3K4me3", True),
    ("H3K4me3/bed", "YD21_k4me3_hg19_fdr_peaks.narrowPeak_rip.csv", False, "H3K4me3", True),
    ("H3K4me3/bed", "YD21_k4me3_hg19_fdr_peaks.narrowPeak_frip.log", False, "H3K4me3", True),

    # Sicer
    ("H3K36me3/bed", "YD16_k36me3_hg19-W200-G1000-FDR1E-6-island.bed", True, "H3K36me3", True),
    ("H3K36me3/bed", "OD16_k36me3_hg19-W200-G1000-FDR1E-6-island.bed", True, "H3K36me3", True),
    ("H3K36me3/bed", "YD18_k36me3_hg19-W200-G1000-FDR1E-6-island.bed_rip.csv", False, "H3K36me3",
     True),

    # All Golden
    ("H3K27ac/bed", "H3K27ac_golden_ODS_weak_consensus.bed", False, "H3K27ac", True),
    ("H3K36me3/bed", "H3K36me3_golden_weak_consensus.bed", False, "H3K36me3", True),

    # # Outliers
    ("H3K27ac/bed", "outliers/YD8_k27ac_hg19_broad_peaks.broadPeak_rip.csv", False, "H3K27ac",
     True),
    ("H3K27ac/bed_all", "YD9_k27ac_hg19_broad_peaks.broadPeak", False, "H3K27ac", True),
    ("H3K27ac/bed_all", "OD9_k27ac_hg19_broad_peaks.broadPeak", False, "H3K27ac", True),
    ("H3K27ac/bed_all", "YD9_k27ac_hg19_broad_peaks.broadPeak", True, "H3K27ac", False),
    ("H3K27ac/bed_all", "OD9_k27ac_hg19_broad_peaks.broadPeak", True, "H3K27ac", False),
    ("H3K27ac/bed_all", "YD8_k27ac_hg19_broad_peaks.broadPeak_rip.csv", False, "H3K27ac", True),
])
def test_collect_golden_peaks(tmp_dir, root, relative_path, selected, mod, exclude_outliers):
    data_root = Path(tmp_dir)

    file = data_root / root / relative_path
    file.parent.mkdir(parents=True)
    file.touch()

    peaks = loi._collect_golden_peaks(data_root, exclude_outliers)

    assert selected == (file in peaks[mod])


@pytest.mark.parametrize("relative_path,selected", [
    ("YD_YD4_H3K27ac_hg19_1.0E-6_peaks.bed", True),
    ("OD_OD7_H3K27ac_hg19_1.0E-6_peaks.bed", True),
    ("YD9_k27ac_hg19_broad_peaks.broadPeak", True),
    ("YD4_k4me3_hg19_fdr_peaks.narrowPeak", True),
    ("YD16_k36me3_hg19-W200-G1000-FDR1E-6-island.bed", True),

    ("OD_OD7_H3K27ac_hg19_1.0E-6_peaks.bed_rip.csv", False),
    ("H3K27ac_zinbra_ODS_weak_consensus.bed", False),
    ("H3K27ac_zinbra_weak_consensus.bed", False),
    ("YD8_k27ac_hg19_broad_peaks.broadPeak_rip.csv", False),
    ("YD21_k4me3_hg19_fdr_peaks.narrowPeak_rip.csv", False),
    ("YD21_k4me3_hg19_fdr_peaks.narrowPeak_frip.log", False),
    ("YD18_k36me3_hg19-W200-G1000-FDR1E-6-island.bed_rip.csv", False),
    ("H3K27ac_golden_ODS_weak_consensus.bed", False),
    ("H3K36me3_golden_weak_consensus.bed", False),
])
def test_collect_peaks_in_folder(tmp_dir, relative_path, selected):
    peaks_root = Path(tmp_dir)

    file = peaks_root / relative_path
    file.touch()

    res = loi._collect_peaks_in_folder(peaks_root)
    assert selected == (file in res)


def generate_test_data_chromhmm(loci_root):
    chromhmm_root = loci_root / "chromhmm"
    chromhmm_root.mkdir()
    (chromhmm_root / "cd14_chromhmm.hg19.10_EnhA2.bed").touch()
    (chromhmm_root / "cd14_chromhmm.hg19.12_ZNF_Rpts.bed").touch()
    (chromhmm_root / "cd14_chromhmm.hg19.9_EnhA1.bed").touch()
    (chromhmm_root / "cd14_chromhmm.hg19.2_TssFlnk.bed").touch()
    (chromhmm_root / "cd14_chromhmm.hg19.1_TssA.bed").touch()
