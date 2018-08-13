import argparse
import os

# run example:
# bash run_python3.sh scripts/python/prepare_methfiles.py
# '-i /scratch/shchukinai/aging/methylation/sorted_data/
# -c /scratch/shchukinai/aging/methylation/cytosines.bed -p methylcall -s minavcov10.sorted
# -o /scratch/shchukinai/aging/methylation/clean_methpipe/data/' methfiles


def convert_line_to_methpipe_format(line):
    id, chr, pos, strand, coverage, freqC, freqT = line.strip().split('\t')
    if strand == 'F':
        strand = '+'
    else:
        strand = '-'
    string = '\t'.join([chr, pos, strand, 'CpG', str(float(freqC) / 100.0), coverage]) + '\n'
    return string


def prepare_methfiles(cytosine_set, folder, prefix, suffix, output_dir):
    # go though input directory and filter all files
    files = [os.path.join(folder, file_name) for file_name in os.listdir(folder)
             if file_name.startswith(prefix)]
    for file_name in files:
        output_name = file_name.split('/')[-1].split('.')[:3]
        output_name.extend([suffix, 'meth'])
        output_name = os.path.join(output_dir, '.'.join(output_name))
        with open(output_name, 'w') as out:
            with open(file_name) as inp:
                is_header = True
                for line in inp:
                    # skip header
                    if is_header:
                        is_header = False
                        continue
                    chr_base = line.split()[0]
                    # skip cytosines that are not in the filtered set
                    if chr_base not in cytosine_set:
                        continue
                    out.write(convert_line_to_methpipe_format(line))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputDirectory', required=True,
                        help='directory with initial sorted files')
    parser.add_argument('-p', '--prefix', help="name prefix for methylation data files",
                        required=True)
    parser.add_argument('-c', '--cytosine', help="BED file with cytosines set", required=True)
    parser.add_argument('-s', '--suffix', help="suffix to add", required=True)
    parser.add_argument('-o', '--outputDirectory', help='output directory')
    args = parser.parse_args()

    directory = args.inputDirectory
    prefix = args.prefix
    set_file = args.cytosine
    output = args.outputDirectory
    suffix = args.suffix

    # read in cytosines to keep
    cytosine_set = set()
    with open(set_file) as inp:
        for line in inp:
            chr, start, end = line.strip().split()
            cytosine_set.add('.'.join([chr, start]))

    prepare_methfiles(cytosine_set, directory, prefix, suffix, output)


if __name__ == "__main__":
    main()
