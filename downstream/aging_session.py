import re
import argparse
from pathlib import Path
from urllib.request import urlopen

# Hardcoded URLs with data
Y20O20_CONSENSUS_PATH = "https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq" \
                        "/Y20O20/peaks/{}/consensus"
Y20O20_BW_PATH = "https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq/Y20O20" \
                 "/bigwigs/{}"
Y20O20_PEAK_PATH = "https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq" \
                   "/Y20O20/peaks/{}"
ENCODE_BW_PATH = "https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq" \
                 "/cd14encode/bigwigs"
ENCODE_PEAK_PATH = "https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq" \
                   "/cd14encode/peaks/{}"
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
TRACK_TEMPLATE = """        <Track altColor="{}" autoScale="false" \
clazz="org.broad.igv.track.FeatureTrack" color="{}" \
colorScale="ContinuousColorScale;0.0;100.0;255,255,255;{}" displayMode="COLLAPSED" \
featureVisibilityWindow="-1" fontSize="10" id="{}" name="{}" renderer="BASIC_FEATURE" \
sortable="false" visible="true" windowFunction="count">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="100.0" \
minimum="0.0" type="LINEAR"/>
        </Track>"""
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
failed_tracks = []


def _cli():
    parser = argparse.ArgumentParser(
        description="For given histone modification generates IGV aging session",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("hist", help="Histone modification")
    parser.add_argument("output", help="Output file path")

    args = parser.parse_args()
    hist = args.hist
    output = args.output

    y20o20_consensuses = search_in_url(Y20O20_CONSENSUS_PATH.format(hist),
                                       '<a href="([^"]*consensus[^"]*)">')
    y20o20_total_consensuses = [y20o20_cons for y20o20_cons in y20o20_consensuses
                                if "OD" not in y20o20_cons and "YD" not in y20o20_cons]
    y20o20_bws = search_in_url(Y20O20_BW_PATH.format(hist), '<a href="([^"]*bw)">')
    y20o20_peaks = search_in_url(Y20O20_PEAK_PATH.format(hist),
                                 '<a href="([^"]*(?:peaks.bed|island.bed|broadPeak))">')
    encode_bws = search_in_url(ENCODE_BW_PATH,
                               '<a href="([^"]*{}[^"]*bw)">'.format(GSM_HIST_MAP[hist]))
    encode_peaks = []
    for tool_path in HIST_TOOL_PATH_MAP[hist]:
        encode_peaks += search_in_url(ENCODE_PEAK_PATH.format(hist) + "/" + tool_path,
                                      '<a href="([^"]*{}[^"]*(?:peaks.bed|island.bed|broadPeak))">'
                                      .format(GSM_HIST_MAP[hist]))

    with open(output, 'w') as f:
        print(HEADER, file=f)

        for y20o20_cons in y20o20_total_consensuses:
            print(RESOURCE_TEMPLATE.format(y20o20_cons), file=f)
        for y20o20_bw in y20o20_bws:
            print(RESOURCE_TEMPLATE.format(y20o20_bw), file=f)
        for y20o20_peak in y20o20_peaks:
            print(RESOURCE_TEMPLATE.format(y20o20_peak), file=f)
        for encode_bw in encode_bws:
            print(RESOURCE_TEMPLATE.format(encode_bw), file=f)
        for encode_peak in encode_peaks:
            print(RESOURCE_TEMPLATE.format(encode_peak), file=f)
        print(RESOURCE_TEMPLATE.format(LABELS_URL.format(hist)), file=f)
        print(RESOURCE_FOOTER, file=f)

        for y20o20_cons in y20o20_total_consensuses:
            color = get_color(hist, y20o20_cons)
            print(TRACK_TEMPLATE.format(color, color, color, y20o20_cons, Path(y20o20_cons).stem),
                  file=f)
        for y20o20_bw in y20o20_bws:
            color = get_color(hist, y20o20_bw)
            print(TRACK_TEMPLATE.format(color, color, color, y20o20_bw,
                                        Path(y20o20_bw).stem.upper()), file=f)
        for y20o20_peak in y20o20_peaks:
            color = get_color(hist, y20o20_peak)
            print(TRACK_TEMPLATE.format(color, color, color, y20o20_peak,
                                        "ZINBRA " + Path(y20o20_peak).stem), file=f)
        for encode_bw in encode_bws:
            color = get_color(hist, encode_bw)
            print(TRACK_TEMPLATE.format(color, color, color, encode_bw, Path(encode_bw).stem),
                  file=f)
        for encode_peak in encode_peaks:
            color = get_color(hist, encode_peak)
            print(TRACK_TEMPLATE.format(color, color, color, encode_peak, Path(encode_peak).stem),
                  file=f)
        print(TRACK_TEMPLATE.format("178,0,0", "178,0,0", "178,0,0", LABELS_URL.format(hist),
                                    Path(LABELS_URL.format(hist)).stem), file=f)

        print(FOOTER, file=f)


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
