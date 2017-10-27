import os
import sys
from pathlib import Path
import pandas as pd
from itertools import chain

# Force matplotlib to not use any Xwindows backend.
import matplotlib
matplotlib.use('Agg')

parent_dir = os.path.dirname(os.path.realpath(__file__))
project_root = os.path.abspath(os.path.join(parent_dir) + "/..")
sys.path.insert(0, project_root)

from reports.bed_metrics import plot_metric_heatmap, bed_metric_table  # noqa
from reports.bed_metrics import heatmap_donor_color_fun as donor_color  # noqa


def process_hist_mod(peaks_root, threads, golden):
    assert golden, "Zinbra paths not supported yet"

    all_peaks_root = peaks_root / "bed_all"

    peaks_paths = list(chain(all_peaks_root.glob("*Peak"),
                             all_peaks_root.glob("*-island.bed")))

    res_prefix = "{}_intersection".format(peaks_root.name)
    df_path = peaks_root / "{}.csv".format(res_prefix)
    if df_path.exists():
        df = pd.DataFrame.from_csv(str(df_path))
    else:
        df = bed_metric_table(peaks_paths, peaks_paths, threads=threads)
        df.to_csv(str(df_path))
        print("Metrics results saved to:", str(df_path))

    def donor_name(name):
        if name.endswith("Peak"):
            tool = "macs2"
        elif name.endswith("island.bed"):
            tool = "sicer"
        else:
            tool = "unknown"
        return name.split('_')[0] + "_" + tool

    result_plot_path = str(peaks_root / "{}.png".format(res_prefix))
    plot_metric_heatmap("Intersection metric: All donors {}".format(
        peaks_root.name
    ),
        df,
        save_to=result_plot_path,
        row_cluster=True, col_cluster=True,
        row_color_fun=donor_color, col_color_fun=donor_color,
        row_label_fun=donor_name, col_label_fun=donor_name
    )
    print("Metrics plot saved to:", str(result_plot_path))


golden_root = Path("/mnt/stripe/bio/experiments/aging/peak_calling")
for mod in golden_root.glob("H*"):
    if mod.is_dir():
        process_hist_mod(Path(mod), threads=30, golden=True)

# if __name__ == "__main__":
