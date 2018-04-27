import argparse
from pathlib import Path

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Script creates pdf report with cumulative consensus plot for selected tools
(all histone modifications).
"""


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("folder", help="Results folder")
    parser.add_argument("output", help="Output pdf path")
    parser.add_argument("tools", nargs='*', help="Tools folders")

    args = parser.parse_args()
    folder = Path(args.folder)
    pdf_path = args.output
    tools = args.tools

    with PdfPages(pdf_path) as pdf:
        for hist_index, hist_mod in enumerate(["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3",
                                               "H3K4me1"]):
            for tool_index, tool in enumerate(tools):
                tool_path = folder / hist_mod / tool

                plt.figure()
                tracks_paths = collect_peaks_in_folder(tool_path)
                if len(tracks_paths) > 0:
                    bar_plot(hist_mod, tool, tracks_paths, pdf, (10, 7), 14)


if __name__ == "__main__":
    # Force matplotlib to not use any Xwindows backend.
    import matplotlib
    matplotlib.use('Agg')

    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    from downstream.loci_of_interest import collect_peaks_in_folder
    from downstream.peak_metrics import bar_plot

    _cli()
