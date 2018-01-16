import argparse
import pandas as pd
from pathlib import Path
from scripts.util import age

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """Script for calculating mean FRIP and peaks count"""


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("peaks", help="Peaks folder")
    args = parser.parse_args()

    tool_outliers = {"H3K27ac": {"zinbra": ["OD16"],  "macs_broad": [], "sicer": []},
                     "H3K27me3": {"zinbra": [],  "macs_broad": ["OD10"], "sicer": []},
                     "H3K36me3": {"zinbra": ["OD20", "OD18"],
                                  "macs_broad": ["OD20", "OD18"],
                                  "sicer": []},
                     "H3K4me3": {"zinbra": ["YD10", "YD15", "YD4", "OD14", "OD6", "OD17"],
                                 "macs_broad": ["YD10", "YD18", "YD5", "OD14", "OD6", "OD17"],
                                 "sicer": []},
                     "H3K4me1": {"zinbra": [],  "macs_broad": ["YD14", "OD7"], "sicer": []}}

    for hist_mod in ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3"]:
        for tool in ["zinbra", "macs_broad", "sicer"]:
            folder_path = Path(args.peaks + "/" + hist_mod + "/" + tool + "/clean")
            rip_files = sorted([str(f) for f in folder_path.glob("*_rip.csv")])

            df = pd.DataFrame(columns=["file", "reads", "peaks", "rip"])
            if len(rip_files) > 0:
                for rip_path in rip_files:
                    data = pd.read_csv(rip_path)
                    rip_file = rip_path.split('/')[-1]
                    reads = float(data["reads"].loc[0])
                    peaks = int(data["peaks"].loc[0])
                    rip = float(data["rip"].loc[0])
                    df.loc[df.size] = (rip_file, reads, peaks, rip)
                df["frip"] = 100.0 * df["rip"] / df["reads"]
                df.index = [age(df.loc[n]["file"]) for n in df.index]

                for donor in tool_outliers[hist_mod][tool]:
                    df = df[df.index != donor]

                sum_frip = 0
                sum_peaks = 0
                mean_frip = df["frip"].mean()
                mean_peaks = df["peaks"].mean()
                for donor in df.index:
                    sum_frip += pow(df["frip"][donor] - mean_frip, 2)
                    sum_peaks += pow(df["peaks"][donor] - mean_peaks, 2)
                print(hist_mod + " " + tool)
                print("FRIP = %0.2f ± %0.2f" % (mean_frip, pow(sum_frip / len(df), 0.5)))
                print("peaks = %0.0f ± %0.0f" % (mean_peaks, pow(sum_peaks / len(df), 0.5)))
                print("")


if __name__ == "__main__":
    _cli()
