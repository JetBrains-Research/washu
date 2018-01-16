import re
from collections import namedtuple

#################################################################
# Age and donors utility code
#################################################################

Age = namedtuple('Age', 'name color prefix')
OLD = Age('O', 'blue', '')
YOUNG = Age('Y', 'red', '')


def is_od_input(c):
    c = re.sub(".*/", "", c)
    input_od = re.match('.*input.*od.*', str(c), flags=re.IGNORECASE) is not None
    od_input = re.match('.*od.*input.*', str(c), flags=re.IGNORECASE) is not None
    return input_od or od_input


def is_yd_input(c):
    c = re.sub(".*/", "", c)
    input_yd = re.match('.*input.*yd.*', str(c), flags=re.IGNORECASE) is not None
    yd_input = re.match('.*yd.*input.*', str(c), flags=re.IGNORECASE) is not None
    return input_yd or yd_input


def is_input(c):
    return is_od_input(c) or is_yd_input(c)


def is_od(c):
    return (re.match('.*od\\d+.*', str(c), flags=re.IGNORECASE) is not None) and not is_input(c)


def is_yd(c):
    return (re.match('.*yd\\d+.*', str(c), flags=re.IGNORECASE) is not None) and not is_input(c)


def is_od_or_yd(c):
    return (re.match('.*[yo]d\\d+.*', str(c), flags=re.IGNORECASE) is not None) and not is_input(c)


def age(n):
    match = re.search('[yo]d\\d+', str(n), flags=re.IGNORECASE)
    return None if not match else match.group(0)


def regions_extension(c):
    return re.match('.*(?:Peak|_peaks\.bed|island\.bed)$', str(c))
