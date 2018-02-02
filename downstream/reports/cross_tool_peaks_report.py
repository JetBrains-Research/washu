import argparse
import datetime
import os
from pathlib import Path

import pandas as pd

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Usage:
    cross_tool_peaks_report.py [input folder] [output folder] [first tool] [second tool]

Script creates pdf reports with peaks statistics:
 1) Bar plot with inclusion and jaccard metrics 
 2) Bar plot with peaks count
"""


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input_folder", help="Data folder")
    parser.add_argument("output_folder", help="Output folder for pdf")
    parser.add_argument("tool1", help="First tool")
    parser.add_argument("tool2", help="Second tool")

    args = parser.parse_args()
    input_folder = Path(args.input_folder)
    output_folder = Path(args.output_folder)
    tool1 = args.tool1
    tool2 = args.tool2

    df_metrics = pd.DataFrame(columns=['hist_mod', 'name', 'type', 'value'])
    df_peaks_count = pd.DataFrame(columns=['hist_mod', 'name', 'type', 'value'])

    today = datetime.datetime.today()
    pdf_path = str(output_folder / (tool1 + "_vs_" + tool2 + "_" + input_folder.parent.name + "_" +
                                    today.strftime("%d.%m.%Y") + ".pdf"))
    for hist_mod in (hist_mod.name for hist_mod in input_folder.glob("*")
                     if os.path.isdir(hist_mod)):
        tracks_paths1 = loi.collect_peaks_in_folder(input_folder / hist_mod / tool1)
        tracks_paths2 = loi.collect_peaks_in_folder(input_folder / hist_mod / tool2)

        if len(tracks_paths1) > 0 and len(tracks_paths2) > 0:
            for track_path1 in tracks_paths1:
                name = track_path1.name.split('_' + hist_mod.lower())[0]
                filtered_tracks_paths2 = [track_path2 for track_path2 in tracks_paths2 if
                                          name in track_path2.name]
                if len(filtered_tracks_paths2) > 0:
                    track_path2 = filtered_tracks_paths2[0]
                    inclusion_first = run_metric_intersection(track_path1, track_path2)[0]
                    inclusion_second = run_metric_intersection(track_path2, track_path1)[0]
                    jaccard = run_metric_jaccard(track_path1, track_path2)[0]
                    peaks_count1 = Bed(str(track_path1)).count()
                    peaks_count2 = Bed(str(track_path2)).count()

                    df_metrics = df_metrics.append({'hist_mod': hist_mod, 'name': name,
                                                    'type': "inclusion_first",
                                                    'value': inclusion_first},
                                                   ignore_index=True)
                    df_metrics = df_metrics.append({'hist_mod': hist_mod, 'name': name,
                                                    'type': "inclusion_second",
                                                    'value': inclusion_second}, ignore_index=True)
                    df_metrics = df_metrics.append({'hist_mod': hist_mod, 'name': name,
                                                    'type': "jaccard", 'value': jaccard},
                                                   ignore_index=True)
                    df_peaks_count = df_peaks_count.append({'hist_mod': hist_mod, 'name': name,
                                                            'type': "first_peaks_count",
                                                            'value': peaks_count1},
                                                           ignore_index=True)
                    df_peaks_count = df_peaks_count.append({'hist_mod': hist_mod, 'name': name,
                                                            'type': "second_peaks_count",
                                                            'value': peaks_count2},
                                                           ignore_index=True)

    with PdfPages(pdf_path) as pdf:
        plot_bar_from_df(df_metrics, pdf)
        plot_bar_from_df(df_peaks_count, pdf)

        desc = pdf.infodict()
        desc['Title'] = 'Report: Combined peaks plots for different callers'
        desc['Author'] = 'JetBrains Research BioLabs'
        desc['Subject'] = 'peaks'
        desc['CreationDate'] = datetime.datetime.today()
        desc['ModDate'] = datetime.datetime.today()


def plot_bar_from_df(df_metrics, pdf):
    plt.figure(figsize=(20, 15))
    g = sns.barplot(data=df_metrics, y="value", x="hist_mod", hue="type", ci="sd")
    g.set_xticklabels(g.get_xticklabels(), rotation=90)
    plt.tight_layout()
    save_plot(pdf)


if __name__ == "__main__":
    # Force matplotlib to not use any Xwindows backend.
    import matplotlib

    matplotlib.use('Agg')

    import seaborn as sns
    import downstream.loci_of_interest as loi
    import matplotlib.pyplot as plt
    from bed.bedtrace import Bed
    from matplotlib.backends.backend_pdf import PdfPages
    from downstream.bed_metrics import save_plot, run_metric_intersection, run_metric_jaccard

    _cli()
