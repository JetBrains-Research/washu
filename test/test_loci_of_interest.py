from pathlib import Path
import pytest
from test.fixtures import test_data, tmp_path
import reports.loci_of_interest as loi


def test_collect_chromhmm(tmp_path):
    loci_root = tmp_path
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


def test_collect_loci(tmp_path):
    loci_root = tmp_path
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
    assert 12 == len(table)
    assert ['chipseq_diff_loci', 'chromhmm', 'else', 'enhancers',
            'golden_consensus', 'golden_median_consensus',
            'regulatory', 'repeats', 'tfs', 'top_level_paths',
            'zinbra_consensus', 'zinbra_median_consensus'
            ] == sorted(str(k) for k in table)
    assert 2 == len(table['chipseq_diff_loci'])
    assert 2 == len(table['top_level_paths'])
    assert 2 == len(table['zinbra_consensus'])
    assert 2 == len(table['zinbra_median_consensus'])
    assert 5 == len(table['chromhmm'])
    assert 27 == sum(len(v) for v in table.values())


@pytest.mark.parametrize("root,relative_path,selected,mod", [
    # Zinbra
    ("peaks/H3K27ac", "YD_YD4_H3K27ac_hg19_1.0E-6_peaks.bed", True, "H3K27ac"),
    ("peaks/H3K27ac", "OD_OD7_H3K27ac_hg19_1.0E-6_peaks.bed", True, "H3K27ac"),
    ("peaks/H3K27ac", "OD_OD7_H3K27ac_hg19_1.0E-6_peaks.bed_rip.csv", False, "H3K27ac"),
    ("peaks/H3K27ac", "H3K27ac_zinbra_ODS_weak_consensus.bed", False, "H3K27ac"),
    ("peaks/H3K27ac", "H3K27ac_zinbra_weak_consensus.bed", False, "H3K27ac"),
    ("peaks/H3K27ac", "H3K27ac_zinbra_weak_consensus.bed", False, "H3K27ac"),
])
def test_collect_zinbra_peaks(tmp_path, root, relative_path, selected, mod):
    data_root = tmp_path
    file = data_root / root / relative_path
    file.parent.mkdir(parents=True)
    file.touch()

    peaks = loi._collect_zinbra_peaks(data_root / "peaks")
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
def test_collect_golden_peaks(tmp_path, root, relative_path, selected, mod, exclude_outliers):
    data_root = tmp_path

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
def test_collect_peaks_in_folder(tmp_path, relative_path, selected):
    peaks_root = tmp_path

    file = peaks_root / relative_path
    file.touch()

    res = loi._collect_peaks_in_folder(peaks_root)
    assert selected == (file in res)


def test_collect_peaks_in_folder_sorted(tmp_path):
    peaks_root = tmp_path

    names = [
        "YD_YD4_H3K27ac_hg19_1.0E-6_peaks.bed",
        "OD_OD7_H3K27ac_hg19_1.0E-6_peaks.bed",
        "YD16_k36me3_hg19-W200-G1000-FDR1E-6-island.bed",
        "OD1_k27ac_hg19_broad_peaks.broadPeak",
        "OD10_k27ac_hg19_broad_peaks.broadPeak",
        "OD4_k4me3_hg19_fdr_peaks.narrowPeak",
        "YD4_k4me3_hg19_fdr_peaks.narrowPeak"
    ]
    files = [peaks_root / name for name in names]
    for f in files:
        f.touch()

    res = loi._collect_peaks_in_folder(peaks_root)
    print([f.name for f in res])
    assert ['OD1_k27ac_hg19_broad_peaks.broadPeak', 'OD4_k4me3_hg19_fdr_peaks.narrowPeak',
            'OD_OD7_H3K27ac_hg19_1.0E-6_peaks.bed', 'OD10_k27ac_hg19_broad_peaks.broadPeak',
            'YD_YD4_H3K27ac_hg19_1.0E-6_peaks.bed', 'YD4_k4me3_hg19_fdr_peaks.narrowPeak',
            'YD16_k36me3_hg19-W200-G1000-FDR1E-6-island.bed'] == [f.name for f in res]


@pytest.mark.parametrize("fname,expected", [
    ("YD_YD4_H3K27ac_hg19_1.0E-6_peaks.bed", ('YD', 4)),
    ("OD_OD7_H3K27ac_hg19_1.0E-6_peaks.bed", ('OD', 7)),
    ("YD9_k27ac_hg19_broad_peaks.broadPeak", ('YD', 9)),
    ("YD4_k4me3_hg19_fdr_peaks.narrowPeak", ('YD', 4)),
    ("YD16_k36me3_hg19-W200-G1000-FDR1E-6-island.bed", ('YD', 16)),
    ("YD4.broadPeak", ('YD', 4)),

    ("H3K27ac_zinbra_ODS_weak_consensus.bed", ('ODS', 'H3K27ac_zinbra_ODS_weak_consensus.bed')),
    ("H3K27ac_zinbra_weak_consensus.bed", ('H3K27ac_zinbra_weak_consensus.bed', 0)),
    ("H3K36me3_golden_weak_consensus.bed", ("H3K36me3_golden_weak_consensus.bed", 0)),
])
def test_donor_order_id(tmp_path, fname, expected):
    assert expected == loi.donor_order_id(tmp_path / fname)


def generate_test_data_chromhmm(loci_root):
    chromhmm_root = loci_root / "chromhmm"
    chromhmm_root.mkdir()
    (chromhmm_root / "cd14_chromhmm.hg19.10_EnhA2.bed").touch()
    (chromhmm_root / "cd14_chromhmm.hg19.12_ZNF_Rpts.bed").touch()
    (chromhmm_root / "cd14_chromhmm.hg19.9_EnhA1.bed").touch()
    (chromhmm_root / "cd14_chromhmm.hg19.2_TssFlnk.bed").touch()
    (chromhmm_root / "cd14_chromhmm.hg19.1_TssA.bed").touch()
