import os
import sys
from pathlib import Path
import pandas as pd

# Force matplotlib to not use any Xwindows backend.
import matplotlib
matplotlib.use('Agg')

parent_dir = os.path.dirname(os.path.realpath(__file__))
project_root = os.path.abspath(os.path.join(parent_dir) + "/..")
sys.path.insert(0, project_root)

from reports.bed_metrics import plot_metric_heatmap, bed_metric_table  # noqa
from reports.bed_metrics import heatmap_donor_color_fun as donor_color

def process_hist_mod(peaks_root, threads):
    all_peaks_root = peaks_root / "bed_all"
    peaks_paths = list(all_peaks_root.glob("*Peak"))

    res_prefix = "{}_intersection".format(peaks_root.name)
    df_path = peaks_root / "{}.csv".format(res_prefix)
    if df_path.exists():
        df = pd.DataFrame.from_csv(str(df_path))
    else:
        df = bed_metric_table(peaks_paths, peaks_paths, threads=threads)
        df.to_csv(str(df_path))
        print("Metrics results saved to:", str(df_path))

    def donor_name(name):
        return name.split('_')[0] + "_golden"

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
        process_hist_mod(Path(mod), threads=30)

# if __name__ == "__main__":
