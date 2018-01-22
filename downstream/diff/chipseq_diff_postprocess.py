import os
import shutil

from downstream.diff.chipseq_diff_plots import DiffProcessor as DiffProcessor

path = "/mnt/stripe/bio/experiments/configs/Y20O20/chip-seq-diff"

output_path = "/mnt/stripe/bio/experiments/configs/Y20O20/chip-seq-diff-locations"


def main():
    for mark in ['H3K27ac', 'H3K27me3', 'H3K4me1', 'H3K4me3', 'H3K36me3']:
        print(mark)
        out = os.path.join(output_path, mark)

        if os.path.exists(out):
            shutil.rmtree(out)
            os.mkdir(out)

        diffs = DiffProcessor(path, out, mark)
        diffs.collect_difference_for_pipeline()


if __name__ == '__main__':
    main()
