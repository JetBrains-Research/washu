from itertools import chain
from pathlib import Path


def collect_loci(loci_root: Path):
    sort_by_fname = lambda p: p.name

    # Chrom HMM
    annotations = {"chromhmm": _collect_chromhmm(loci_root)}

    # Other Loci
    top_level_paths = []
    for f in loci_root.iterdir():
        if f.is_dir():
            if (f.name != "chromhmm"):
                annotations[f.name] = sorted(f.glob('*.bed'), key=sort_by_fname)
        else:
            # add top level files
            top_level_paths.append(f)

    # All annotations:
    annotations[None] = sorted(chain(top_level_paths, *annotations.values()),
                               key=sort_by_fname)

    # Default annotations: top level + selected folders
    default_paths = list(top_level_paths)
    for key in ["enhancers", "tfs", "regulatory", "repeats",
                "golden_consensus", "zinbra_consensus"]:
        default_paths.extend(annotations[key])
    annotations["default"] = sorted(default_paths, key=sort_by_fname)

    return annotations


def chromhmm_state_descr(fname):
    """
    :param fname: State file name
    :return: State human readable descriptor
    """
    chunks = fname.split(".")
    if len(chunks) > 2:
        state = chunks[2]
        if state in _CHROMHMM_ST_MAP:
            return "{} ({})".format(_CHROMHMM_ST_MAP.get(chunks[2]), chunks[2])
    return fname


_CHROMHMM_ST_MAP = {
    "1_TssA": "Active TSS",
    "2_TssFlnk": "Flanking TSS",
    "3_TssFlnkU": "Flanking TSS Upstream",
    "4_TssFlnkD": "Flanking TSS Downstream",
    "5_Tx": "Strong transcription",
    "6_TxWk": "Weak transcription",
    "7_EnhG1": "Genic enhancer1",
    "8_EnhG2": "Genic enhancer2",
    "9_EnhA1": "Active Enhancer 1",
    "10_EnhA2": "Active Enhancer 2",
    "11_EnhWk": "Weak Enhancer",
    "12_ZNF_Rpts": "ZNF genes & repeats",
    "13_Het": "Heterochromatin",
    "14_TssBiv": "Bivalent/Poised TSS",
    "15_EnhBiv": "Bivalent Enhancer",
    "16_ReprPC": "Repressed PolyComb",
    "17_ReprPCWk": "Weak Repressed PolyComb",
    "18_Quies": "Quiescent/Low",
}


def _collect_chromhmm(loci_root):
    return sorted((loci_root / "chromhmm").glob('*.bed'),
                  key=lambda p: int(p.name.split(".")[2].split("_")[0]))


def _collect_peaks(peaks_roots, outliers=False):
    result = {}
    for peaks_root in [x for x in peaks_roots.iterdir() if x.is_dir() and x.name.startswith("H")]:
        print("Peaks:", peaks_root)

        peaks = list(chain(peaks_root.glob("**/*_peaks.bed"),
                           peaks_root.glob("**/bed/*-island.bed"),
                           peaks_root.glob("**/bed/*.*Peak")
                           # peaks_root.glob("**/*consensus*.bed")  # Ignore
                          ))
        for p in peaks:
            assert "outlier" not in str(p)
        # e.g.
        # * OD_OD14_H3K27ac_hg19_1.0E-6_peaks.bed
        # * OD8_k27ac_hg19_broad_peaks.broadPeak
        # * Ignore consensus: e.g. zinbra_weak_consensus.bed
        #TODO: peaks.sort(key=donor_order_id)
        print(len(peaks))
        print(*[str(p) for p in peaks], sep="\n")
        result[peaks_root.name] = peaks
    return result


if __name__ == "__main__":
    # data_root = Path("/Volumes/BigData/bio")
    data_root = Path("/mnt/stripe/bio")

    loci_root = data_root / "raw-data/aging/loci_of_interest"
    golden_peaks_root = data_root / "experiments/aging/peak_calling"
    zinbra_peaks_root = data_root / "experiments/configs/Y20O20/peaks"

    # Full:
    #zinbra_peaks_root = data_root / "experiments/configs/Y20O20_full/peaks"

    signal_root = data_root / "experiments/signal"
