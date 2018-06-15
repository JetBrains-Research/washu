import argparse
import pandas as pd
from pathlib import Path
from downstream.aging import donor

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """Script for calculating mean FRIP and peaks count"""
failed_tracks_path = "/mnt/stripe/bio/experiments/aging/Y20O20.failed_tracks.csv"
failed_tracks_df = pd.read_csv(failed_tracks_path, delimiter="\t", skiprows=1, index_col="donor")


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("peaks", help="Peaks folder")
    args = parser.parse_args()

    tool_outliers = {"H3K27ac": {"span": [],  "macs_broad": [], "sicer": []},
                     "H3K27me3": {"span": ["OD10"],  "macs_broad": ["OD10"], "sicer": []},
                     "H3K36me3": {"span": ["OD18"],
                                  "macs_broad": ["OD18"],
                                  "sicer": []},
                     "H3K4me3": {"span": ["YD18", "YD5", "OD17"],
                                 "macs_broad": ["YD18", "YD5", "OD17"],
                                 "sicer": []},
                     "H3K4me1": {"span": [],  "macs_broad": ["YD14", "OD7"], "sicer": []}}

    for hist_mod in ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3"]:
        for tool in ["span", "macs_broad", "sicer"]:
            print(hist_mod, tool)
            folder_path = Path(args.peaks + "/" + hist_mod + "/" + tool)
            rip_files = sorted([str(f) for f in folder_path.glob("*_rip.csv")])

            df = pd.DataFrame(columns=["file", "reads", "peaks", "length", "rip"])
            if len(rip_files) > 0:
                for rip_path in rip_files:
                    data = pd.read_csv(rip_path).fillna(0)
                    rip_file = rip_path.split('/')[-1]
                    reads = float(data["reads"].loc[0])
                    peaks = int(data["peaks"].loc[0])
                    length = int(data["length"].loc[0])
                    rip = float(data["rip"].loc[0])
                    df.loc[df.size] = (rip_file, reads, peaks, length, rip)
                df["frip"] = 100.0 * df["rip"] / df["reads"]
                df.index = [donor(df.loc[n]["file"]) for n in df.index]

                for donor in failed_tracks_df.loc[:, hist_mod].index:
                    if failed_tracks_df.loc[:, hist_mod][donor] == 1:
                        df = df[df.index != donor]

                for donor in tool_outliers[hist_mod][tool]:
                    df = df[df.index != donor]

                sum_frip = 0
                sum_peaks = 0
                sum_length = 0
                mean_frip = df["frip"].mean()
                mean_peaks = df["peaks"].mean()
                mean_length = df["length"].mean()
                for donor in df.index:
                    sum_frip += pow(df["frip"][donor] - mean_frip, 2)
                    sum_peaks += pow(df["peaks"][donor] - mean_peaks, 2)
                    sum_length += pow(df["length"][donor] - mean_length, 2)
                print("FRIP = %0.2f ± %0.2f" % (mean_frip, pow(sum_frip / len(df), 0.5)))
                print("peaks = %0.0f ± %0.0f" % (mean_peaks, pow(sum_peaks / len(df), 0.5)))
                print("length = %0.0f ± %0.0f" % (mean_length, pow(sum_length / len(df), 0.5)))
                print("")


if __name__ == "__main__":
    _cli()
