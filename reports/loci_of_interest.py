from itertools import chain
from pathlib import Path


def collect_loci(loci_root: Path):
    sort_by_fname = lambda p: p.name

    # Chrom HMM
    annotations = {}

    # Other Loci
    top_level_paths = []
    for f in loci_root.iterdir():
        if f.is_dir():
            if (f.name != "chromhmm"):
                annotations[f.name] = sorted(f.glob('**/*.bed'), key=sort_by_fname)
            else:
                annotations[f.name] = sorted(
                    f.glob('*.bed'),
                    key=lambda p: int(p.name.split(".")[2].split("_")[0])
                )
        elif f.suffix == ".bed":
            # add top level files
            top_level_paths.append(f)

    # All annotations:
    annotations[None] = sorted(chain(top_level_paths, *annotations.values()),
                               key=sort_by_fname)

    # Default annotations: top level + selected folders (~ FIXED, doesn't depend on to peak
    # calling)
    default_paths = list(top_level_paths)
    for key in ["enhancers", "tfs", "regulatory", "repeats", 'chromhmm']:
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


def _collect_peaks_in_folder(peaks_root):
    return sorted(chain(peaks_root.glob("*_peaks.bed"),
                        peaks_root.glob("*-island.bed"),
                        peaks_root.glob("*.*Peak")))


def _collect_zinbra_peaks(zinbra_peaks_root):
    zinbra_peaks = {}
    for folder in (zinbra_peaks_root / "peaks").iterdir():
        if folder.is_dir() and folder.name.startswith("H"):
            zinbra_peaks[folder.name] = _collect_peaks_in_folder(folder)

    return zinbra_peaks


def _collect_golden_peaks(golden_peaks_root, exclude_outliers):
    golden_peaks = {}
    for folder in golden_peaks_root.iterdir():
        if folder.is_dir() and folder.name.startswith("H"):
            sub_folder = "bed" if exclude_outliers else "bed_all"
            golden_peaks[folder.name] = _collect_peaks_in_folder(folder / sub_folder)
    return golden_peaks
