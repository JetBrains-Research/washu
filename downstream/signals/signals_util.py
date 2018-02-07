import re
from pathlib import Path


def extract_normalization(path: Path):
    folder = path.parent.name
    file_wo_ext = path.stem

    if file_wo_ext.startswith(folder):
        return file_wo_ext[len(folder) + 1:]

    return file_wo_ext


def extract_datatype(path: Path):
    matches = re.match(".*/(H[a-z0-9]+)/.*", str(path), re.IGNORECASE)

    if matches:
        return matches.group(1)

    if "/meth/" in str(path):
        return "meth"

    return "N/A"
