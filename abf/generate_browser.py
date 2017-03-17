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


# TRACKS example
# <!--K27ac YD-->
#                   {name:                 'K27ac YD5',
#                    bwgURI:               '/k27ac_10vs10_bams_bws/YD5_R1_hg19.bw',
#                    style: [{type: 'default', style: {glyph: 'HISTOGRAM', BGCOLOR: 'rgb(219, 41, 35)', HEIGHT: 30, id: 'style1'}}],
#                    noDownsample:         true},

BrowserRecord = namedtuple('BrowserRecord', ['name', 'browser_type', 'dirs'])


def process_record(output, record):
    """Creates corresponding html with name and content and returns <li> element for landing page"""
    print("Processing", record)
    if record.browser_type != 'bw':
        print("Unknown type", record.browser_type)
        # TODO[shpynov] Add folders browsing
        return None

    tracks = []
    for dir in record.dirs.split(','):
        print("Processing", dir)
        for f in glob.glob('{}/*.bw'.format(dir)):
            print("BW", f)
            # TODO[shpynov] added colors
            tracks.append("""
                   {name:                 '$NAME',
                    bwgURI:               '/$FILE',
                    style: [{type: 'default', style: {glyph: 'HISTOGRAM', BGCOLOR: 'rgb(219, 41, 35)', HEIGHT: 30, id: 'style1'}}],
                    noDownsample:         true},
""".replace('$NAME', record.name).replace('$FILE', f))
    with open('{}/data.html'.format(path), 'r') as file:
        data_html = file.read()
    file = '{}/{}.html'.format(output, record.name)
    with open(file, 'w') as file:
        file.write(data_html.
                   replace('$DATA_NAME', record.name).
                   replace('//$DATA_HERE', '\n'.join(tracks)))
    print('Created', file)
    return '<li><a href="{0}/{1}.html">{1}</a></li>'.format(output, record.name)


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
generate_browser.py is a script to generate index.html and other html files for ageing web server
USAGE: generate_browser.py --browsers name1 type1 folders11,...,folder1n1 name2 type2 folder21,..,folder2n2 output_folder
This will generate index html with 2 links name1 and name2 containing embedded biodallance browsers for containing bw files.
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
