import datetime
import argparse
import tempfile
import numpy as np
import pandas as pd

from pathlib import Path

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Usage: consensuses_range_report.py [peaks summary] [output pdf path] [top peaks count (optional)]
[histone modification] [consensus percent] [tools list]

Script creates pdf report with ChIP-seq consensus statistics:
 1) peak consensus venn diagram
 2) peak consensus bar plot
"""


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("peaks_summary", help="Peaks summary")
    parser.add_argument("output", help="Output pdf path")
    parser.add_argument('-p', '--peaks_count', help="Top peaks count (optional)",
                        type=int, default=0)
    parser.add_argument('-t', '--threads', help="Threads number for parallel processing",
                        type=int, default=6)
    parser.add_argument("hist_mod", help="Histone modification")
    parser.add_argument("percent", help="Consensus percent")
    parser.add_argument("tools", nargs='*', help="Tools folders")

    args = parser.parse_args()
    peaks_summary = args.peaks_summary
    pdf_path = args.output
    hist_mod = args.hist_mod
    percent = int(args.percent)
    peaks_count = args.peaks_count
    threads_num = args.threads
    tools = args.tools

    tmp_dir = Path(tempfile.gettempdir())
    dfd = pd.read_csv(peaks_summary, sep='\t', comment='#')
    dfd = dfd[['modification', 'tool', 'procedure', 'status', 'file']]

    with PdfPages(pdf_path) as pdf:
        for tool_index, tool in enumerate(tools):
            curr_dfd = dfd.loc[np.logical_and(dfd['modification'] == hist_mod,
                                              dfd['tool'] == tool)]
            if len(curr_dfd) > 0:
                for index, procedure in enumerate(["tuned", "default"]):
                    tracks_paths = curr_dfd.loc[
                        np.logical_and(curr_dfd['procedure'] == procedure,
                                       curr_dfd['status'] == 'ok')]['file']

                    filtered_paths = []
                    if peaks_count != 0:
                        for track_path in tracks_paths:
                            tmp_path = tmp_dir / "{}_{}.bed".format(Path(track_path).stem,
                                                                    peaks_count)
                            with open(str(tmp_path), 'w') as f:
                                run((["sort", "-k9nr", str(track_path)],
                                     ["head", "-n", str(peaks_count)],
                                     ["sort", "-k1,1", "-k2,2n", "-k3,3n"]),
                                    stdout=f)
                                filtered_paths.append(tmp_path)
                    else:
                        filtered_paths = tracks_paths

                    od_paths_map = {donor(track_path): track_path for track_path in filtered_paths
                                    if regions_extension(track_path) and is_od(track_path)}
                    yd_paths_map = {donor(track_path): track_path for track_path in filtered_paths
                                    if regions_extension(track_path) and is_yd(track_path)}

                    # Code for different consensuses investigation
                    od_consensus_bed, yd_consensus_bed, yd_od_int_bed = \
                        calc_consensus_file(list(od_paths_map.values()),
                                            list(yd_paths_map.values()),
                                            percent=percent)
                    venn_consensus(od_consensus_bed, yd_consensus_bed, percent, pdf,
                                   title_prefix=hist_mod + " " + tool + " " + procedure)
                    bar_consensus(od_paths_map, yd_paths_map, od_consensus_bed,
                                  yd_consensus_bed, yd_od_int_bed, threads_num, pdf,
                                  percent=percent,
                                  title_prefix=hist_mod + " " + tool + " " + procedure)

        desc = pdf.infodict()
        desc['Title'] = 'Report: Consensus plots for data investigation'
        desc['Author'] = 'JetBrains Research BioLabs'
        desc['Subject'] = 'consensus'
        desc['CreationDate'] = datetime.datetime.today()
        desc['ModDate'] = datetime.datetime.today()


if __name__ == "__main__":
    # Force matplotlib to not use any Xwindows backend.
    import matplotlib
    matplotlib.use('Agg')

    from matplotlib.backends.backend_pdf import PdfPages
    from downstream.aging import regions_extension, donor, is_od, is_yd
    from bed.bedtrace import run
    from downstream.peak_metrics import calc_consensus_file, venn_consensus, bar_consensus

    _cli()
