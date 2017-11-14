import subprocess
from pipeline_utils import run_bash, PROJECT_ROOT_PATH

import pytest
from test.fixtures import test_data


@pytest.mark.parametrize("help_arg", [True, False])
def test_consensus_help(capfd, help_arg):
    if help_arg:
        args = ["-s", "foo", "-h"]
    else:
        args = []

    with pytest.raises(subprocess.CalledProcessError) as e:
        run_bash("bed/consensus.sh", *args)

    out, _err = capfd.readouterr()
    assert "returned non-zero exit status 1" in str(e.value)
    assert """Calculate consensus for peaks in selected folder

Usage: consensus.sh [OPTIONS] folder_path

Options:
  -p number       Consensus percent should be taken as number (cannot be used with -c)
  -c number       Count of tracks for consensus should be taken as number
                          (cannot be used with -p)
  -h|--help       Show help

With no arguments - show help""" in out
