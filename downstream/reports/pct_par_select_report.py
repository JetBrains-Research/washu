import sys
import os
import argparse
from pathlib import Path

def donor_order_id(path):
    name = path.stem.split('_')[1]
    return (name, int(name[2:]))

def calc_overlap(folder_path):
    peaks_paths = sorted(folder_path.glob("*peaks.bed"),
                         key=donor_order_id)
    return bm.bed_metric_table(peaks_paths, peaks_paths, threads=30)


def strategy_folders(folder_path):
    return os.walk(folder_path)[0][1]

def _cli():
    parser = argparse.ArgumentParser(description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("folder")
    parser.add_argument("output")
    args = parser.parse_args()
    folder_path = Path(args.folder)
    output_path = Path(args.output)
    folders = strategy_folders(folder_path)
    print("Calculating overlap")
    df = pd.DataFrame(columns=['strategy','a','b','overlap'])
    for path in folders:
        overlap = calc_overlap(path)
        for a in overlap.index:
            for b in overlap.columns:
                df.loc[len(df)] = (path.name, a, b, df[a, b])
    df.to_csv(output)


if __name__ == "__main__":
    import downstream.bed_metrics as bm

    _cli()