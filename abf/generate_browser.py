#!/usr/bin/env python
import argparse
import glob
import shutil
from collections import namedtuple
import os

path = os.path.dirname(os.path.realpath(__file__))

__author__ = 'Oleg Shpynov'


class WritableDirectory(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if not os.path.isdir(values):
            raise argparse.ArgumentTypeError("{0} is not a valid path".format(values))
        if os.access(values, os.R_OK + os.W_OK):
            setattr(namespace, self.dest, values)
        else:
            raise argparse.ArgumentTypeError("{0} is not a writable directory".format(values))


def mod_color(file):
    if 'k27ac' in file:
        return 255, 0, 0
    if 'k4me1' in file:
        return 231, 45, 56
    return 0, 0, 178


BrowserRecord = namedtuple('BrowserRecord', ['name', 'browser_type', 'dirs'])


def process_record(output, record):
    """Creates corresponding html with name and content and returns <li> element for landing page"""
    print("Processing", record)
    if record.browser_type == 'bw':
        tracks = []
        for dir in record.dirs.split(','):
            print("Processing", dir)
            for f in glob.glob('{}/*.bw'.format(dir)):
                print("BW", f)
                r, g, b = mod_color(f)
                # OD are shown 20% darker
                if 'OD' in f:
                    r, g, b = int(0.8 * r), int(0.8 * g), int(0.8 * b)
                rgb = "{},{},{}".format(r, g, b)

                tracks.append("""
{name:                 '$NAME',
bwgURI:               '/$FILE',
style: [{type: 'default', style: {glyph: 'HISTOGRAM', BGCOLOR: 'rgb($RGB)', HEIGHT: 30, id: 'style1'}}],
noDownsample:         true},
""".replace('$NAME', os.path.basename(f)).replace('$FILE', f).replace('$RGB', rgb))
        with open('{}/data.html'.format(path), 'r') as file:
            data_html = file.read()
        file = '{}/{}.html'.format(output, record.name)
        with open(file, 'w') as file:
            file.write(data_html.
                       replace('$DATA_NAME', record.name).
                       replace('//$DATA_HERE', '\n'.join(tracks)))
        print('Created', file)
        return '<li><a href="{0}/{1}.html">{1}</a></li>'.format(output, record.name)
    elif record.browser_type == 'dir':
        ds = record.dirs.split(',')
        if len(ds) > 1:
            raise argparse.ArgumentTypeError('Type dir should have single dir')
        return '<li><a href="{0}">{1}</a></li>'.format(ds[0], record.name)
    else:
        print("Unknown type", record.browser_type)
        return None


def process(output, records):
    # Read in the file
    with open('{}/index.html'.format(path), 'r') as file:
        index = file.read()
    browsers = '\n'.join(filter(lambda x: x is not None, [process_record(output, r) for r in records]))
    index_file = '{}/index.html'.format(output)
    with open(index_file, 'w') as file:
        file.write(index.replace('$CONTENT_HERE', browsers))
    print('Created', index_file)
    shutil.copy('{}/style.css'.format(path), output)


def main():
    parser = argparse.ArgumentParser(description='''
generate_browser.py is a script to generate index.html and other html files for ageing web server.

USAGE: generate_browser.py --output out --browsers name1 type1 folders11,...,folder1n1 name2 type2 folder21,..,folder2n2
Supported types: bw and dir.

EXAMPLE:
python ~/work/washu/abf/generate_browser.py --output . --browsers k27ac bw k27ac_10vs10_bams_bws \
    k27ac_signal bw k27ac_10vs10_signal \
    k27ac_peaks dir k27ac_10vs10_bams_macs_broad_0.01 \
    k4me1 bw k4me1_10vs10_reseq_bams_bws \
    k4me1_signal bw k4me1_10vs10_reseq_signal \
    k4me1_peaks dir k4me1_10vs10_reseq_bams_macs_broad_0.01
''')
    parser.add_argument('--output', required=True, action=WritableDirectory, type=str,
                        help='Path to directory with data to run pipeline')
    parser.add_argument('--browsers', required=True, metavar='name browser_type dirs', type=str, nargs='+',
                        help='browsers to create in format: name type dirs')
    args = parser.parse_args()
    if len(args.browsers) % 3 != 0:
        raise argparse.ArgumentTypeError()

    records = []
    it = iter(args.browsers)
    for name, browser_type, dirs in zip(it, it, it):
        records.append(BrowserRecord(name, browser_type, dirs))
    print(records)
    process(args.output, records)


if __name__ == "__main__":
    main()
