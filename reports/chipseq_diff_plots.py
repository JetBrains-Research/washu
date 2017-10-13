import os
import shutil
import tempfile
import pandas as pd
import glob
from matplotlib_venn import venn2, venn3


# Some utility functions to analyse ChIP-seq difference

def count_lines(file):
    with open(file) as f:
        lines = f.readlines()
        return len(lines)


class DiffProcessor:
    def __init__(self, input, output, mark):
        self.input_path = os.path.join(input, mark)
        self.output = output
        self.mark = mark
        self.diff_counts = []
        self.bed_files_produced = []

    def process_diff_bind(self):
        diff_bind_path = os.path.join(self.input_path, "{}_diffbind".format(self.mark))
        for file in glob.glob(os.path.join(diff_bind_path, '*_difference.csv')):
            base_name = os.path.basename(file)
            df = pd.read_csv(file)
            self.diff_counts.append(("diffbind_" + base_name, df.shape[0]))

    def add_bed_file(self, base_name, file):
        self.diff_counts.append((base_name, count_lines(file)))
        shutil.copyfile(file, os.path.join(self.output, base_name))
        self.bed_files_produced.append(base_name)

    def process_zinbra(self):
        zinbra_base_path = "/mnt/stripe/bio/experiments/configs/Y20O20/enrichment"
        for file in glob.glob(os.path.join(zinbra_base_path, "diff_OD_YD_{}_zinbra*".format(self.mark))):
            base_name = os.path.basename(file)
            if "cond" in base_name:
                continue

            self.diff_counts.append((base_name, count_lines(file)))

            os.system("cut -f1-3 {} > {}".format(file, os.path.join(self.output, base_name)))

            self.bed_files_produced.append(base_name)

    def process_macs_bg_diff(self):
        folder_name = "{}_macs_bdgdiff".format(self.mark)
        macs_bdgdiff = os.path.join(self.input_path, folder_name)
        size = 0
        for file in glob.glob(os.path.join(macs_bdgdiff, "{}_*_cond*.bed".format(self.mark))):
            size += count_lines(file) - 1
        self.diff_counts.append((folder_name, size))

    def process_macs_pooled(self):
        folder_name = "{}_macs_pooled".format(self.mark)
        macs_pooled = os.path.join(self.input_path, folder_name)

        size = 0
        result_base_file_name = "{}.bed".format(folder_name)
        with open(os.path.join(self.output, result_base_file_name), "w") as out:
            for file in glob.glob(os.path.join(macs_pooled, "{}_*_cond*.bed".format(self.mark))):
                size += count_lines(file)

                with open(file) as f:
                    out.writelines(f.readlines())

        self.bed_files_produced.append(result_base_file_name)

        self.diff_counts.append((folder_name, size))

    def process_macs_pooled_Y_vs_O(self):
        folder_name = "{}_macs_pooled_Y_vs_O".format(self.mark)
        macs_pooled_Y_vs_O = os.path.join(self.input_path, folder_name)

        size = 0
        result_base_file_name = "{}.bed".format(folder_name)
        with open(os.path.join(self.output, result_base_file_name), "w") as out:
            for file in glob.glob(os.path.join(macs_pooled_Y_vs_O, "{}_*.broadPeak".format(self.mark))):
                size += count_lines(file)

                with open(file) as f:
                    out.writelines(f.readlines())

        self.diff_counts.append((folder_name, size))
        self.bed_files_produced.append(result_base_file_name)

    def process_diffreps(self):
        folder_name = os.path.join("/mnt/stripe/bio/experiments/configs/Y20O20/diffreps", self.mark)

        self.add_bed_file("diffreps_{}.bed".format(self.mark), os.path.join(folder_name, "hotspot.bed"))

    def process_diffreps_not_uniqe(self):
        folder_name = os.path.join("/mnt/stripe/kurbatsky/configs/Y20O20/diffreps", self.mark)

        self.add_bed_file("diffreps_not_uniqe_{}.bed".format(self.mark), os.path.join(folder_name, "hotspot.bed"))

    def compare_difference(self):

        self.process_diff_bind()

        self.process_zinbra()

        self.process_macs_bg_diff()

        self.process_macs_pooled()

        self.process_macs_pooled_Y_vs_O()

        self.process_diffreps()

        self.process_diffreps_not_uniqe()

    def count_intersections(self):
        temp_dir = tempfile.mkdtemp(suffix=".tmp")
        print(temp_dir)

        union_sh = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                                 '../bed/union.sh'))

        result = {}

        for p1 in self.bed_files_produced:
            for p2 in self.bed_files_produced:
                union_file = os.path.join(temp_dir, "counts.bed")
                command = "bash {} {} {} >{}".format(
                    union_sh,
                    os.path.join(self.output, p1),
                    os.path.join(self.output, p2),
                    union_file)
                os.system(command)

                intersection = 0.0

                with open(union_file) as f:
                    for l in f.readlines():
                        if "|" in l:
                            intersection += 1

                result[(p1, p2)] = intersection

        shutil.rmtree(temp_dir)
        return result

    def plot_venn2(self, f1, f2):
        temp_dir = tempfile.mkdtemp(suffix=".tmp")

        union_sh = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                                 '../bed/union.sh'))

        result = [0.0] * 3

        union_file = os.path.join(temp_dir, "counts.bed")
        command = "bash {} {} {} >{}".format(
            union_sh,
            os.path.join(self.output, f1),
            os.path.join(self.output, f2),
            union_file)
        print(command)
        os.system(command)

        with open(union_file) as f:
            for l in f.readlines():
                index = 0
                if f1 in l: index += 1
                if f2 in l: index += 2
                result[index - 1] += 1

        shutil.rmtree(temp_dir)

        venn2(subsets=result, set_labels=(f1, f2))

    def plot_venn3(self, f1, f2, f3):
        temp_dir = tempfile.mkdtemp(suffix=".tmp")

        union_sh = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                                 '../bed/union.sh'))

        result = [0.0] * 7

        union_file = os.path.join(temp_dir, "counts.bed")
        command = "bash {} {} {} {} >{}".format(
            union_sh,
            os.path.join(self.output, f1),
            os.path.join(self.output, f2),
            os.path.join(self.output, f3),
            union_file)
        print(command)
        os.system(command)

        with open(union_file) as f:
            for l in f.readlines():
                index = 0
                if f1 in l: index += 1
                if f2 in l: index += 2
                if f3 in l: index += 4
                result[index - 1] += 1


        shutil.rmtree(temp_dir)

        venn3(subsets=result, set_labels=(f1, f2, f3))

