import re
from collections import namedtuple
from scripts.util import is_input

Group = namedtuple('Group', 'name color')
OLD = Group('O', 'blue')
YOUNG = Group('Y', 'red')
UNKNOWN = Group('', 'black')


def is_od(c):
    return not is_input(c) and group(c) == OLD


def is_yd(c):
    return not is_input(c) and group(c) == YOUNG


def is_od_or_yd(c):
    return is_od(c) or is_yd(c)


def group(c):
    if re.match('.*od\\d+.*', str(c), flags=re.IGNORECASE):
        return OLD
    if re.match('.*yd\\d+.*', str(c), flags=re.IGNORECASE):
        return YOUNG
    return UNKNOWN


def donor(c):
    name = re.sub('.*/', '', str(c))
    match = re.search('[yo]d\\d+', name, flags=re.IGNORECASE)
    return name if not match else match.group(0)


def regions_extension(c):
    return re.match('.*(?:(?:broad|narrow)Peak|_peaks(?:_\d+)?\.bed|island(?:_\d+)?\.bed)$', str(c))
