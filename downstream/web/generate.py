import os
import re
import shutil
from pathlib import Path

# Hardcoded URLs with data
from downstream.sessions.aging_session import get_color, search_in_url

Y20O20_CONSENSUS_PATH = "https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq" \
                        "/Y20O20/peaks/{}/zinbra/consensus"
Y20O20_BW_PATH = "https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq/Y20O20" \
                 "/bigwigs/{}"
Y20O20_BB_PATH = "https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq" \
                 "/Y20O20/peaks/{}/{}/bigBed"
ENCODE_BW_PATH = "https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq" \
                 "/cd14encode/bigwigs"
ENCODE_BB_PATH = "https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq" \
                 "/cd14encode/peaks/{}/{}/bigBed"
LABELS_URL = "https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq/Y20O20" \
             "/labels/{}_labels.bed"

GSM_HIST_MAP = {
    'H3K27ac': 'GSM1102782',
    'H3K27me3': 'GSM1102785',
    'H3K36me3': 'GSM1102788',
    'H3K4me1': 'GSM1102793',
    'H3K4me3': 'GSM1102797'
}
HIST_TOOL_PATH_MAP = {
    'H3K27ac': {'macs_broad'},
    'H3K27me3': {'macs_broad', 'sicer'},
    'H3K36me3': {'macs_broad', 'sicer'},
    'H3K4me1': {'macs_broad'},
    'H3K4me3': {'macs_broad'}
}
HIST_COLOR_MAP = {
    "H3K27ac": (255, 0, 0),
    "H3K27me3": (153, 0, 255),
    "H3K4me1": (255, 153, 0),
    "H3K4me3": (51, 204, 51),
    "H3K36me3": (0, 0, 204)
}
TOOL_COLOR_MAP = {
    "broad": (55, 126, 184),
    "island": (228, 26, 28),
    "peaks": (152, 78, 163),
    "zinbra": (152, 78, 163)
}

HREF_MASK = '<a href="([^"]*.{})">'

failed_tracks = []

FOLDER = os.path.dirname(__file__)
OUT_FOLDER = FOLDER + '/out'


def create_hist_page(hist, hist_page):
    print('Creating modification page {}'.format(hist_page))
    tracks = []

    y20o20_consensuses = search_in_url(Y20O20_CONSENSUS_PATH.format(hist),
                                       '<a href="([^"]*consensus[^"]*)">')
    y20o20_total_consensuses = [y20o20_cons for y20o20_cons in y20o20_consensuses
                                if "OD" not in y20o20_cons and "YD" not in y20o20_cons]
    y20o20_bws = search_in_url(Y20O20_BW_PATH.format(hist), HREF_MASK.format("bw"))
    y20o20_zinbra_peaks = search_in_url(Y20O20_BB_PATH.format(hist, 'zinbra'),
                                        HREF_MASK.format("bb"))
    encode_bws = search_in_url(ENCODE_BW_PATH,
                               '<a href="([^"]*{}[^"]*.bw)">'.format(GSM_HIST_MAP[hist]))
    encode_peaks = []
    for tool_path in HIST_TOOL_PATH_MAP[hist]:
        encode_peaks += search_in_url(ENCODE_BB_PATH.format(hist, tool_path),
                                      '<a href="([^"]*{}[^"]*.bb)">'.format(GSM_HIST_MAP[hist]))

        insensitive_hist = re.compile(re.escape(hist), re.IGNORECASE)

        tracks.extend(print_tracks(hist, y20o20_total_consensuses))
        for i in range(0, len(y20o20_bws)):
            tracks.extend(print_tracks(hist, [y20o20_bws[i]],
                                       name_processor=lambda x: insensitive_hist.sub(hist, x.upper())))
            tracks.extend(print_tracks(hist, [y20o20_zinbra_peaks[i]],
                                       name_processor=lambda x: "ZINBRA " + x))

        tracks.extend(print_tracks(hist, encode_bws))
        tracks.extend(print_tracks(hist, encode_peaks,
                                   name_processor=lambda x: ("MACS " if "broad" in x else "SICER ") + x))
        tracks.append(format_track(hist, LABELS_URL.format(hist)))

    print('Loading data template')
    with open(FOLDER + '/_data.html', 'r') as file:
        data_template_html = file.read()
    with open(OUT_FOLDER + '/' + hist_page, 'w') as file:
        file.write(data_template_html
                   .replace('//@TRACKS@', ',\n'.join(tracks))
                   .replace('@MODIFICATION@', hist.upper()))


def print_tracks(hist, paths, name_processor=lambda x: x):
    return [format_track(hist, path, name_processor) for path in paths]


def format_track(hist, uri, name_processor=lambda x: x):
    return """{
    name: '@NAME@',
    bwgURI: '@URI@',
    style: [{type: 'default', style: {glyph: 'HISTOGRAM', BGCOLOR: 'rgb(@RGB@)', HEIGHT: 30, id: 'style1'}}],
    noDownsample: false}""". \
        replace('@NAME@', name_processor(Path(uri).stem)).replace('@URI@', uri).replace('@RGB@', get_color(hist, uri))


def _cli():

    if os.path.exists(OUT_FOLDER):
        shutil.rmtree(OUT_FOLDER)
    os.mkdir(OUT_FOLDER)
    print('Generating site structure in {}'.format(OUT_FOLDER))

    print('Copying resources')
    for file in os.listdir(FOLDER):
        if re.match('.*\.(html|css|png)', file) and not re.match('template\\.html', file):
            shutil.copy(file, OUT_FOLDER)

    print('Loading template')
    with open(FOLDER + '/template.html', 'r') as file:
        template_html = file.read()

    print('Creating index.html')
    with open(OUT_FOLDER + '/index.html', 'w') as file:
        file.write(template_html.
                   replace('@TITLE@', 'Aging project').
                   replace('@SCRIPTS@', '').
                   replace('@CONTENT@', '_index.html'))

    print('Creating paper.html')
    with open(OUT_FOLDER + '/paper.html', 'w') as file:
        file.write(template_html.
                   replace('@TITLE@', 'Paper').
                   replace('@SCRIPTS@', '').
                   replace('@CONTENT@', '_paper.html'))

    print('Creating software.html')
    with open(OUT_FOLDER + '/software.html', 'w') as file:
        file.write(template_html.
                   replace('@TITLE@', 'Software').
                   replace('@SCRIPTS@', '').
                   replace('@CONTENT@', '_software.html'))

    print('Creating explore pages')
    for hist in GSM_HIST_MAP.keys():
        hist_page = 'data_{}.html'.format(hist).lower()
        create_hist_page(hist, hist_page)
        # biodallance should be included within main html page
        with open(OUT_FOLDER + '/{}.html'.format(hist).lower(), 'w') as file:
            file.write(template_html.
                       replace('@TITLE@', 'hist').
                       replace('@SCRIPTS@',
                               """<script type="text/javascript" src="//www.biodalliance.org/release-0.13/dalliance-compiled.js"></script>""").
                       replace('@CONTENT@', hist_page))


if __name__ == "__main__":
    _cli()
