import re
import argparse
from pathlib import Path
from urllib.request import urlopen

# Hardcoded URLs with data
ENCODE_BW_PATH = "https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq" \
                 "/cd14encode/bigwigs"
ENCODE_PEAK_PATH = "https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq" \
                   "/cd14encode/peaks/{}"
ENCODE_BB_PATH = "https://www.encodeproject.org/files/{}/@@download/{}.bigBed"
LABELS_URL = "https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq/Y20O20" \
             "/labels/{}_labels.bed"
#'Encode GSM1102782 macs narrow'
ENCODE_HIST_MAP = {
    'H3K27ac': {'GSM': 'GSM1102782', 'ENC': 'ENCFF439NLA'},
    'H3K27me3': {'GSM': 'GSM1102785', 'ENC': 'ENCFF575VMI'},
    'H3K36me3': {'GSM': 'GSM1102788', 'ENC': 'ENCFF003KSH'},
    'H3K4me1': {'GSM': 'GSM1102793', 'ENC': 'ENCFF158WJC'},
    'H3K4me3': {'GSM': 'GSM1102797', 'ENC': 'ENCFF317WLK'}
}
HIST_TOOL_PATH_MAP = {
    'H3K27ac': {'macs_broad/1.0E-4'},
    'H3K27me3': {'macs_broad/1.0E-4', 'sicer/1.0E-6_600'},
    'H3K36me3': {'macs_broad/1.0E-4', 'sicer/1.0E-6_600'},
    'H3K4me1': {'macs_broad/1.0E-4'},
    'H3K4me3': {'macs_broad/1.0E-4'}
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
HEADER = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="hg19" hasGeneTrack="true" hasSequenceTrack="true" locus="All" \
path="/home/user/aging/${HIST}_igv_session.xml" version="8">
    <Resources>"""
RESOURCE_TEMPLATE = "        <Resource path=\"{}\"/>"
RESOURCE_FOOTER = """    </Resources>
    <Panel height="591" name="DataPanel" width="1836">"""
IGV_TRACK_TEMPLATE = """        <Track altColor="{}" autoScale="true" \
clazz="org.broad.igv.track.FeatureTrack" color="{}" \
colorScale="ContinuousColorScale;0.0;100.0;255,255,255;{}" displayMode="COLLAPSED" \
featureVisibilityWindow="-1" fontSize="10" id="{}" name="{}" renderer="BASIC_FEATURE" \
sortable="false" visible="true" windowFunction="count">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="100.0" \
minimum="0.0" type="LINEAR"/>
        </Track>"""
UCSC_TRACK_TEMPLATE = """track name={} {} visibility={} maxHeightPixels=50:100:11 \
windowingFunction=maximum smoothingWIndow=off {} {}"""
FOOTER = """    </Panel>
    <Panel height="367" name="FeaturePanel" width="1836">
        <Track altColor="0,0,178" autoScale="false" color="0,0,178" displayMode="COLLAPSED" \
featureVisibilityWindow="-1" fontSize="10" id="Reference sequence" name="Reference sequence" \
sortable="false" visible="true"/>
        <Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" \
color="0,0,178" colorScale="ContinuousColorScale;0.0;423.0;255,255,255;0,0,178" \
displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="35" id="hg19_genes" \
name="RefSeq Genes" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="423.0" \
minimum="0.0" type="LINEAR"/>
        </Track>
    </Panel>
    <PanelLayout dividerFractions="0.6145077720207254"/>
    <HiddenAttributes>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
        <Attribute name="NAME"/>
    </HiddenAttributes>
</Session>"""
IGV_BROWSER = "igv"
UCSC_BROWSER = "ucsc"
failed_tracks = []


def _cli():
    parser = argparse.ArgumentParser(
        description="For given histone modification generates IGV aging session",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("hist", help="Histone modification")
    parser.add_argument("output", help="Output file path")
    parser.add_argument("browser", help="Browser name (igv or ucsc)")

    args = parser.parse_args()
    hist = args.hist
    output = args.output
    browser = args.browser

    encode_bws = search_in_url(ENCODE_BW_PATH,
                               '<a href="([^"]*{}[^"]*bw)">'.format(ENCODE_HIST_MAP[hist]['GSM']))
    encode_default_peaks = []
    for tool_path in HIST_TOOL_PATH_MAP[hist]:
        encode_default_peaks += search_in_url(
            ENCODE_PEAK_PATH.format(hist) + "/" + tool_path,
            '<a href="([^"]*{}[^"]*(?:peaks.bed|island.bed|broadPeak))">'
            .format(ENCODE_HIST_MAP[hist]['GSM']))
    encode_tuned_peaks = search_in_url(ENCODE_PEAK_PATH.format(hist),
        '<a href="([^"]*{}[^"]*(?:peaks.bed|island.bed|broadPeak))">'.format(
            ENCODE_HIST_MAP[hist]['GSM']))

    with open(output, 'w') as f:
        big_bed_path = ENCODE_BB_PATH.format(ENCODE_HIST_MAP[hist]['ENC'],
                                             ENCODE_HIST_MAP[hist]['ENC'])

        if browser == IGV_BROWSER:
            print(HEADER, file=f)

            for path in encode_bws + encode_default_peaks + encode_tuned_peaks + \
                        [big_bed_path] + [LABELS_URL.format(hist)]:
                print(RESOURCE_TEMPLATE.format(path), file=f)
            print(RESOURCE_FOOTER, file=f)

        print_tracks(hist, browser, encode_bws, "type=bigWig", f, visibility="full")
        print_tracks(hist, browser, encode_default_peaks, "", f, visibility="dense")
        print(format_track(browser, ','.join(str(x) for x in TOOL_COLOR_MAP["broad"]),
                           big_bed_path, "type=bigBed", visibility="dense",
                           name_processor=lambda _:
                           "Encode {} macs narrow".format(ENCODE_HIST_MAP[hist]['GSM'])), file=f)
        print_tracks(hist, browser, encode_tuned_peaks, "type=broadPeak", f, visibility="dense")
        print(format_track(browser, "", LABELS_URL.format(hist), "", visibility="dense"), file=f)

        if browser == IGV_BROWSER:
            print(FOOTER, file=f)


def print_tracks(hist, browser, paths, track_type, file, visibility="",
                 name_processor=lambda x: x):
    for path in paths:
        color = get_color(hist, path)
        print(format_track(browser, color, path, track_type, visibility=visibility,
                           name_processor=name_processor), file=file)


def format_track(browser, color, path, track_type, visibility="", name_processor=lambda x: x):
    if browser == IGV_BROWSER:
        return IGV_TRACK_TEMPLATE.format(color, color, color, path, name_processor(Path(path).stem))
    if browser == UCSC_BROWSER:
        track_type = "type=broadPeak" if "broadPeak" in path else track_type
        return UCSC_TRACK_TEMPLATE.format(name_processor(Path(path).stem),
                                          "" if urlopen(path).getheader("Content-Length") == '0'
                                          else track_type, visibility,
                                          "itemRgb=On" if "labels" in path else "color=" + color,
                                          ("bigDataUrl=" if "big" in track_type else "\n") + path)


def get_color(hist, filename):
    if "failed" in filename.lower() or re.match('.*[YO]D\\d+.*', filename, flags=re.IGNORECASE) \
            and re.findall('[YO]D\\d+', filename, flags=re.IGNORECASE)[0].lower() in failed_tracks:
        failed_tracks.append(re.findall('[YO]D\\d+', filename, flags=re.IGNORECASE)[0].lower())
        return "192,192,192"
    color = HIST_COLOR_MAP[hist]
    for tool in TOOL_COLOR_MAP.keys():
        if tool in filename:
            color = TOOL_COLOR_MAP[tool]
            break
    if "od" in filename.lower():
        return ','.join(str(int(x * 0.7)) for x in color)
    else:
        return ','.join(str(x) for x in color)


def search_in_url(url, regexp):
    html = str(urlopen(url).read())
    file_names = re.findall(regexp, html)
    file_urls = [url + "/" + file_name for file_name in file_names]
    return file_urls


if __name__ == "__main__":
    _cli()
