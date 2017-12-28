from itertools import chain
from pathlib import Path
import sys


def collect_loci(loci_root: Path):
    sort_by_fname = lambda p: p.name

    # Chrom HMM
    annotations = {}

    # Other Loci
    top_level_paths = []
    for f in loci_root.iterdir():
        if f.is_dir():
            if ("chromhmm" in f.name):
                annotations[f.name] = sorted(
                    f.glob('*.bed'),
                    key=lambda p: int(p.name.split(".")[2].split("_")[0])
                )
            else:
                annotations[f.name] = sorted(f.glob('**/*.bed'), key=sort_by_fname)
        elif f.suffix == ".bed":
            # add top level files
            top_level_paths.append(f)

    annotations["top_level_paths"] = sorted(top_level_paths, key=sort_by_fname)

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
            return "{} ({})".format(chunks[2], _CHROMHMM_ST_MAP.get(chunks[2]))
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


def donor_order_id(path):
    chunks = path.stem.split('_')
    cands = list(filter(lambda s: len(s) > 2 and (s.startswith("OD") or s.startswith("YD")),
                        chunks))
    if len(cands) > 0:
        donor_id = cands[0]
        if donor_id[2] != "S":
            return (donor_id[:2], int(donor_id[2:]))
        else:
            return (donor_id[:3], path.name)
    return (path.name, 0)


def label_converter_shorten_loci(name):
    if "chromhmm" in name:
        idx = name.index("chromhmm")
        suffix = name[:idx]
        return suffix + chromhmm_state_descr(name[idx:])

    name = name.replace(".bed", "")
    name = name.replace("median_consensus", "mcs")
    name = name.replace("weak_consensus", "wcs")
    name = name.replace("without", "w/o")
    return name


def collect_peaks_in_folder(peaks_root):
    return sorted(chain(peaks_root.glob("*_peaks.bed"),
                        peaks_root.glob("*-island.bed"),
                        peaks_root.glob("*.*Peak")),
                  key=donor_order_id)


def _collect_zinbra_peaks(zinbra_peaks_root):
    zinbra_peaks = {}
    for folder in zinbra_peaks_root.iterdir():
        if folder.is_dir() and folder.name.startswith("H"):
            zinbra_peaks[folder.name] = collect_peaks_in_folder(folder)

    return zinbra_peaks


def _collect_golden_peaks(golden_peaks_root, exclude_outliers):
    golden_peaks = {}
    for folder in golden_peaks_root.iterdir():
        if folder.is_dir() and folder.name.startswith("H"):
            sub_folder = "bed" if exclude_outliers else "bed_all"
            # TODO: temporary solution while directories layout not fixed
            if (folder / sub_folder).exists():
                # legacy
                golden_peaks[folder.name] = collect_peaks_in_folder(folder / sub_folder)
            else:
                # new version
                golden_peaks[folder.name] = collect_peaks_in_folder(folder)
    return golden_peaks
