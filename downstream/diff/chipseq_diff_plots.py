import os
import shutil
import tempfile
import pandas as pd
import glob
import subprocess

import sys
from matplotlib_venn import venn2, venn3


# Some utility functions to analyse ChIP-seq difference

def count_lines(file):
    with open(file) as f:
        lines = f.readlines()
        return len(lines)


def count_if_exists(bed):
    if os.path.exists(bed):
        return count_lines(bed)
    else:
        return 0


class ChangeCollector:
    def __init__(self, input, output, mark, change_type):
        assert change_type in ["young", "old", "both"]

        self.input_path = os.path.join(input, mark)

        if not os.path.exists(output):
            os.makedirs(output)

        self.output = output

        self.mark = mark
        self.change_type = change_type
        self.change_counts = []

        self.change_files_produced = []

    def get_file_name(self, result_name):
        return result_name + "_" + self.change_type + ".bed"

    def process_diff_bind(self, name):
        diff_bind_path = os.path.join(self.input_path, name)
        for file in glob.glob(os.path.join(diff_bind_path, '*_difference.csv')):
            df = pd.read_csv(file)
            result_name = "{}_{}".format(name, os.path.splitext(os.path.basename(file))[0])
            if len(df) == 0.0:
                self.change_counts.append((result_name, 0))
                continue

            if self.change_type == "young":
                df = df[df["Fold"] < 0.0]

            if self.change_type == "old":
                df = df[df["Fold"] > 0.0]

            self.change_counts.append((result_name, df.shape[0]))

            result_file_name = self.get_file_name(result_name)
            with open(os.path.join(self.output, result_file_name), "w") as f:
                for i, row in df.iterrows():
                    f.write("{}\t{}\t{}\n".format(row["seqnames"], row["start"], row["end"]))

            self.change_files_produced.append(result_file_name)

    def add_bed_file(self, result_name, file):
        self.change_counts.append((result_name, count_lines(file)))
        result_file_name = self.get_file_name(result_name)
        shutil.copyfile(file, os.path.join(self.output, result_file_name))
        self.change_files_produced.append(result_file_name)

    def process_zinbra(self, input):
        base_path = "/mnt/stripe/bio/experiments/configs/Y20O20"
        zinbra_name = "zinbra_input" if input else "zinbra"
        zinbra_base_path = "{}/chip-seq-diff/{}/{}".format(base_path, self.mark, zinbra_name)
        pattern = os.path.join(zinbra_base_path, "diff_OD_YD_{}_zinbra*".format(self.mark))
        for file in glob.glob(pattern):
            base_name = os.path.basename(file)
            if "cond" in base_name:
                continue

            if self.change_type == "young":
                file = os.path.splitext(file)[0] + "_cond2.bed"

            if self.change_type == "old":
                file = os.path.splitext(file)[0] + "_cond1.bed"

            result_name = os.path.splitext(base_name)[0]
            self.change_counts.append((result_name, count_if_exists(file)))

            if os.path.exists(file):
                result_file_name = self.get_file_name(result_name)

                os.system("cut -f1-3 {} > {}".format(file, os.path.join(self.output,
                                                                        result_file_name)))

                self.change_files_produced.append(result_file_name)

    def process_macs_bg_diff(self):
        folder_name = "macs_bdgdiff"
        macs_bdgdiff = os.path.join(self.input_path, folder_name)
        size = 0
        for file in glob.glob(os.path.join(macs_bdgdiff, "{}_*_cond*.bed".format(self.mark))):
            size += count_lines(file) - 1
        self.change_counts.append((folder_name + "_" + self.mark, size))

    def process_diffreps(self, folder_name):
        folder = os.path.join(self.input_path, folder_name)

        if self.change_type == "both":
            name = "enriched.bed"
        elif self.change_type == "young":
            name = "enriched_down.bed"
        elif self.change_type == "old":
            name = "enriched_up.bed"
        else:
            raise ValueError("Wrong change_type: {}".format(self.change_type))

        self.add_bed_file("{}_{}".format(folder_name, self.mark),
                          os.path.join(folder, name))

        if self.change_type == "both":
            self.add_bed_file("{}_{}_hotspot".format(folder_name, self.mark),
                              os.path.join(folder, "hotspot.bed"))
        else:
            self.change_counts.append(("{}_{}_hotspot".format(folder_name, self.mark), 0))

    def process_chip_diff(self):
        folder = os.path.join(self.input_path, "chipdiff")

        file1 = os.path.join(folder, "{}_3_cond1.bed".format(self.mark))
        file2 = os.path.join(folder, "{}_3_cond2.bed".format(self.mark))
        if self.change_type == "both":
            size = count_lines(file1) + count_lines(file2)
        elif self.change_type == "young":
            size = count_lines(file1)
        elif self.change_type == "old":
            size = count_lines(file2)
        else:
            raise ValueError()

        self.change_counts.append(("chipdiff_{}".format(self.mark), size))

    def collect_difference(self):
        self.process_diff_bind("diff_bind_zinbra")
        self.process_diff_bind("diff_bind_cons_zinbra")

        self.process_chip_diff()

        self.process_zinbra(False)
        self.process_zinbra(True)

        self.process_macs_bg_diff()

        self.process_diffreps("diffReps")
        self.process_diffreps("diffReps_broad")
        self.process_diffreps("diffReps_broad_input")
        self.process_diffreps("diffReps_input")

    def count_intersections(self):
        temp_dir = tempfile.mkdtemp(suffix=".tmp")
        print(temp_dir)

        union_sh = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                                 '../../bed/union.sh'))

        result = {}

        for p1 in self.change_files_produced:
            for p2 in self.change_files_produced:
                union_file = os.path.join(temp_dir, "counts.bed")
                command = "bash {} {} {} >{}".format(
                    union_sh,
                    os.path.join(self.output, p1),
                    os.path.join(self.output, p2),
                    union_file)
                run_result = subprocess.run(command, shell=True, stderr=subprocess.PIPE)
                if run_result.returncode != 0:
                    print(run_result)
                    print(run_result.stderr, file=sys.stderr)
                    raise Exception("command failed: {}".format(command))

                intersection = 0.0

                with open(union_file) as f:
                    for l in f.readlines():
                        if "|" in l:
                            intersection += 1

                result[(p1, p2)] = intersection

        shutil.rmtree(temp_dir)
        return result

    def collect_difference_for_pipeline(self):
        self.process_macs_pooled()

        self.process_diffreps("diffReps_broad")
        self.process_diffreps("diffReps_broad_input")


class DiffProcessor:
    def __init__(self, input, output, mark):

        self.young = ChangeCollector(input, os.path.join(output, "young"), mark, "young")
        self.old = ChangeCollector(input, os.path.join(output, "old"), mark, "old")
        self.both = ChangeCollector(input, os.path.join(output, "both"), mark, "both")

        self.input_path = os.path.join(input, mark)

        if not os.path.exists(output):
            os.makedirs(output)

        self.output = output

        self.mark = mark

    def get_counts(self):
        result = []
        for (file, count) in self.both.change_counts:
            result.append([file, count])

        young_dict = dict(self.young.change_counts)

        for r in result:
            name = r[0]
            r.append(young_dict[name])

        old_dict = dict(self.old.change_counts)

        for r in result:
            name = r[0]
            r.append(old_dict[name])

        return result

    def collect_difference(self):
        self.young.collect_difference()
        self.old.collect_difference()
        self.both.collect_difference()

    def collect_difference_for_pipeline(self):
        self.young.collect_difference_for_pipeline()
        self.old.collect_difference_for_pipeline()
        self.both.collect_difference_for_pipeline()

    def count_intersections(self):
        return self.both.count_intersections()

    def get_bed_files_produced(self):
        return self.both.change_files_produced

    def plot_venn2(self, f1, f2):
        temp_dir = tempfile.mkdtemp(suffix=".tmp")

        union_sh = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                                 '../../bed/union.sh'))

        result = [0.0] * 3

        union_file = os.path.join(temp_dir, "counts.bed")
        command = "bash {} {} {} >{}".format(
            union_sh,
            os.path.join(self.both.output, f1),
            os.path.join(self.both.output, f2),
            union_file)
        print(command)
        status = os.system(command)
        print(status)

        with open(union_file) as f:
            for l in f.readlines():
                index = 0
                if f1 in l:
                    index += 1
                if f2 in l:
                    index += 2
                result[index - 1] += 1

        shutil.rmtree(temp_dir)

        venn2(subsets=result, set_labels=(f1, f2))

    def plot_venn3(self, f1, f2, f3):
        temp_dir = tempfile.mkdtemp(suffix=".tmp")

        union_sh = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                                 '../../bed/union.sh'))

        result = [0.0] * 7

        union_file = os.path.join(temp_dir, "counts.bed")
        command = "bash {} {} {} {} >{}".format(
            union_sh,
            os.path.join(self.both.output, f1),
            os.path.join(self.both.output, f2),
            os.path.join(self.both.output, f3),
            union_file)
        print(command)
        os.system(command)

        with open(union_file) as f:
            for l in f.readlines():
                index = 0
                if f1 in l:
                    index += 1
                if f2 in l:
                    index += 2
                if f3 in l:
                    index += 4
                result[index - 1] += 1

        shutil.rmtree(temp_dir)

        venn3(subsets=result, set_labels=(f1, f2, f3))
