from pipeline_utils import *


def expand_region(start, end, min_length, chr_length):
    delta = (min_length - (end - start + 1)) / 2
    start = round(max(0, start - delta))
    end = round(min(chr_length, end + delta))
    return start, end


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--sorted_cytosines', required=True,
                        help='set of cytosines (has to be sorted)')
    parser.add_argument('-l', '--readLength', required=True,
                        help='size of the read to use when merging')
    parser.add_argument('-o', '--outputFile', help='output file name', required=True)
    parser.add_argument('-e', '--expand', help='whether to expand regions to the read length',
                        action='store_true')
    parser.add_argument('-i', '--path_to_indexes', action=WritableDirectory, type=str,
                        help='Path to indexes')
    parser.add_argument('-g', '--genome', type=str, help='Genome')
    args = parser.parse_args()

    all_cytosines = args.sorted_cytosines
    read_length = int(args.readLength)
    output_file = args.outputFile
    to_expand = args.expand
    genome = args.genome
    indexes = os.path.join(args.path_to_indexes, genome)

    chrom_sizes_file = os.path.join(indexes, genome + ".chrom.sizes")

    # read chromosome sizes to not exceed length
    chrom_sizes = {}
    with open(chrom_sizes_file) as inp:
        for line in inp:
            chr, length = line.strip().split()
            chrom_sizes[chr] = int(length)

    # compile covered regions: merge cytosines that are closer than read_length
    results = []
    with open(all_cytosines) as inp:
        is_initialized = False
        for line in inp:
            chr, start, end = line.strip().split()
            start, end = int(start), int(end)
            if not is_initialized:
                cur_chr, cur_start, cur_end = chr, start, end
                is_initialized = True
                continue
            if chr == cur_chr and start - read_length <= cur_end:
                cur_end = end
            else:
                # expand small regions: region can't be smaller than read_length if flag is set
                if to_expand and cur_end - cur_start + 1 < read_length:
                    cur_start, cur_end = expand_region(cur_start, cur_end, read_length,
                                                       chrom_sizes[cur_chr])
                results.append((cur_chr, cur_start, cur_end))
                cur_chr, cur_start, cur_end = chr, start, end

    # write results to file
    with open(output_file, 'w') as out:
        for (chr, start, end) in results:
            out.write('\t'.join([chr, str(start), str(end)]) + '\n')


if __name__ == "__main__":
    main()
