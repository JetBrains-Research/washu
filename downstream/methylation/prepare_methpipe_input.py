import argparse
import os
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputDataDirectory', help='directory with .meth files',
                        required=True)
    parser.add_argument('-n', '--name', help="name to use for table and design", required=True)
    parser.add_argument('-o', '--outputDirectory', help='output directory')
    args = parser.parse_args()

    data_directory = args.inputDataDirectory
    name_prefix = args.name
    output_directory = args.outputDirectory

    # use Methpipe to prepare table
    table_name = os.path.join(output_directory, name_prefix + '_proportion_table.txt')
    methfiles = [os.path.join(data_directory, file_name) for file_name
                 in os.listdir(data_directory)]
    # subprocess.run('module load gcc-4.7.2', shell=True)
    subprocess.run('merge-methcounts -t ' + ' '.join(methfiles) + ' > {}'.format(table_name),
                   shell=True)

    # construct design table using proportion table header
    with open(table_name) as inp:
        prop_header = inp.readline().strip().split()
    design_name = os.path.join(output_directory, name_prefix + '_design_matrix_complex.txt')
    with open(design_name, 'w') as out:
        out.write('\t'.join(['base', 'case', 'batch']) + '\n')
        base = '1'
        for methfile in prop_header:
            case = '0'
            batch = '0'
            id = methfile.split('.')[2]
            if id.startswith('OD'):
                case = '1'
            if int(id[2:]) > 10:
                batch = '1'
            out.write('\t'.join([methfile, base, case, batch]) + '\n')


if __name__ == "__main__":
    main()
