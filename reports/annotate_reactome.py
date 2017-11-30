import argparse
import requests
import contextlib
import logging
import time
import sys
import re

import lxml.html
from lxml import etree
from pathlib import Path
import pandas as pd


def _lxml_try(*, sleep_s):
    def deco(fun):
        @contextlib.wraps(fun)
        def inner(*args, sleep_s=sleep_s, **kw):
            for i in range(10):
                try:
                    return fun(*args, **kw)
                except (requests.exceptions.RequestException, etree.LxmlError) as ex:
                    msg = "{}:{} An exception of type {} occurred. " \
                          "Args:\n{!r}"
                    logging.warning(msg.format(id, i, type(ex).__name__,
                                               ex.args))

                    # Wait and retry
                    if sleep_s:
                        time.sleep(sleep_s)
            raise etree.LxmlError("Cannot fetch data for {}".format(id))

        return inner

    return deco


@_lxml_try(sleep_s=0.01)
def fetch_title(s):
    reactome_id = extract_reactome_id(s)

    if not reactome_id:
        return "N/A"
    else:
        url = "http://www.reactome.org/content/detail/" + reactome_id
        t = lxml.html.parse(url)
        return t.find(".//title").text.replace("Reactome | ", "")


def extract_reactome_id(s):
    match = re.match(".*R-HSA-(\d+).*", str(s), flags=re.IGNORECASE)
    if match:
        reactome_id = match.group(1)
    else:
        reactome_id = None
    return reactome_id


def _cli():
    parser = argparse.ArgumentParser(
        description="Annotates table containing http://www.reactome.org pathways "
                    "ids (R-HSA-nnnnn) with pathways titles",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('path', metavar="PATH",
                        help="Table path")

    parser.add_argument('col', metavar="INT_OR_STR",
                        help="Index (>=1) or name of column containing Reactome pathway id. If "
                             "column is and index column w/o name, use ''")

    parser.add_argument('-F', dest="sep", metavar="FS",
                        help="Field separator, auto-detect if not specified")

    parser.add_argument('-o', dest="output_path", metavar="PATH",
                        help="Output path")
    args = parser.parse_args()

    path = Path(args.path)
    sep = args.sep
    col = args.col
    output_path = args.output_path

    try:
        col_i = int(col)
        if col_i <= 0:
            print("Column index should be >= 1, but was:", col, file=sys.stderr)
            sys.exit(-1)
        col_idx = col_i - 1
    except ValueError:
        col_idx = None

    if not path.exists():
        print("File not exists:", str(path), file=sys.stderr)
        sys.exit(1)

    if not path.is_file():
        print("File isn't a regular file:", str(path), file=sys.stderr)
        sys.exit(1)

    if path.stat().st_size == 0:
        print("File is empty: ", str(path), file=sys.stderr)
        sys.exit(0)

    index_col = 0 if col_idx is None and col == "" else None
    df = pd.read_csv(str(path), sep=sep, header='infer' if col_idx is None else None,
                     index_col=index_col)
    if col_idx is None:
        print("Table header:", ", ".join("'{}'".format(c) for c in df.columns), file=sys.stderr)
        if index_col is not None:
            reactome_ids = df.index
        else:
            reactome_ids = df.loc[:, col]
    else:
        reactome_ids = df.iloc[:, col_idx]

    titles = []
    print("Fetch titles for:", len(reactome_ids), "records", file=sys.stderr)
    for i, s in enumerate(reactome_ids, 1):
        titles.append(fetch_title(s))
        if i % 100 == 0:
            print('\n{} of {}'.format(i, len(reactome_ids)), file=sys.stderr)
        else:
            print('.', end="", file=sys.stderr)
            sys.stderr.flush()
    print(file=sys.stderr)

    df['reactome_titles'] = titles
    output_kw = dict(index=(index_col is not None),
                     header=True if col_idx is None else False)
    if output_path:
        print("Save to: ", output_path, file=sys.stderr)
        if sep:
            output_kw['sep'] = sep
        df.to_csv(output_path, **output_kw)
    else:
        # pd.set_option('display.large_repr', 'truncate')
        # pd.set_option('display.max_columns', 0)
        pd.set_option('max_colwidth', 800)
        print(df.to_string(**output_kw))


if __name__ == "__main__":
    _cli()
