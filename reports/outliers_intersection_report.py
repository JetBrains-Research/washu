import os
import sys
from pathlib import Path
import pandas as pd

parent_dir = os.path.dirname(os.path.realpath(__file__))
project_root = os.path.abspath(os.path.join(parent_dir) + "/..")
sys.path.insert(0, project_root)

from reports.bed_metrics import plot_metric_heatmap, bed_metric_table  # noqa


def process_hist_mod(peaks_root, threads):
    all_peaks_root = peaks_root / "bed_all"
    peaks_paths = list(all_peaks_root.glob("*Peak"))

    df_path = all_peaks_root / "intersection.csv"
    if df_path.exists():
        df = pd.DataFrame.from_csv(str(df_path))
    else:
        df = bed_metric_table(peaks_paths, peaks_paths, threads=threads)
        df.to_csv(str(df_path))

    print(df.head())
    plot_metric_heatmap("Intersection metric: All donors {}".format(
        peaks_root.name
    ), save_to=str(all_peaks_root / "intersection.png"))


h3k4me3 = Path("/mnt/stripe/bio/experiments/aging/peak_calling/H3K4me3")
process_hist_mod(h3k4me3, threads=30)

# if __name__ == "__main__":
