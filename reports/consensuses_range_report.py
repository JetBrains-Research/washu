import os
import sys
import tempfile
import datetime
from pathlib import Path

# Force matplotlib to not use any Xwindows backend.
import matplotlib

matplotlib.use('Agg')

parent_dir = os.path.dirname(os.path.realpath(__file__))
project_root = os.path.abspath(os.path.join(parent_dir) + "/..")
sys.path.insert(0, project_root)

from matplotlib.backends.backend_pdf import PdfPages  # nopep8
from scripts.util import regions_extension, age, is_od, is_yd, is_od_or_yd  # nopep8
from bed.bedtrace import Bed, run  # nopep8
from reports.peak_metrics import calc_consensus, venn_consensus, bar_consensus  # nopep8

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Usage: peak_metrics.py [peaks folder] [output pdf path] [top peaks count (optional)]

Script creates pdf report with ChIP-seq consensus statistics:
 1) peak consensus (15%, 20%, 33%, 50%, 66%, 80%, 85%, 100%) venn diagram
 2) peak consensus (15%, 20%, 33%, 50%, 66%, 80%, 85%, 100%) bar plot
"""


def main():
    args = sys.argv

    if len(args) < 2:
        print(help_data)
        sys.exit(1)

    folder_path = Path(args[1])
    paths = sorted([str(f) for f in folder_path.iterdir() if regions_extension(f.name)])
    tmp_dir = Path(tempfile.gettempdir())
    filtered_paths = []

    if len(args) == 4:
        for bed_path in paths:
            tmp_path = tmp_dir / "{}_{}.bed".format(Path(bed_path).stem, args[3])
            with open(str(tmp_path), 'w') as f:
                run((["sort", "-k9nr", str(bed_path)], ["head", "-n", args[3]]), stdout=f)
                filtered_paths.append(tmp_path.name)
    else:
        filtered_paths = paths

    tracks_paths = sorted({path for path in filtered_paths if is_od_or_yd(path)})
    od_paths_map = {age(track_path): Bed(track_path) for track_path in tracks_paths
                    if regions_extension(track_path) and is_od(track_path)}
    yd_paths_map = {age(track_path): Bed(track_path) for track_path in tracks_paths
                    if regions_extension(track_path) and is_yd(track_path)}

    with PdfPages(args[2]) as pdf:
        # Code for different consensuses investigation
        for scale in [1, 1.1667, 1.25, 1.5, 2, 3, 5, 7]:
            od_consensus_bed, yd_consensus_bed, yd_od_int_bed = \
                calc_consensus(od_paths_map, yd_paths_map, scale)
            venn_consensus(od_consensus_bed, yd_consensus_bed, scale, pdf)
            bar_consensus(od_paths_map, yd_paths_map, od_consensus_bed, yd_consensus_bed,
                          yd_od_int_bed, num_of_threads, pdf)

        desc = pdf.infodict()
        desc['Title'] = 'Report: Consensus plots for data investigation'
        desc['Author'] = 'JetBrains Research BioLabs'
        desc['Subject'] = 'consensus'
        desc['CreationDate'] = datetime.datetime.today()
        desc['ModDate'] = datetime.datetime.today()


if __name__ == "__main__":
    num_of_threads = 30
    main()
