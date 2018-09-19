__author__ = 'Roman Chernyatchik'
__email__ = 'roman.chernyatchik@jetbrains.com'

import os
import pytest
from pathlib import Path
from bed.bedtrace import _cleanup


@pytest.fixture
def test_data():
    def inner(relative_path):
        td = os.path.dirname(os.path.abspath(__file__)) + '/testdata'
        return os.path.join(td, relative_path)

    return inner


@pytest.fixture(autouse=True)
def bedtrace_cleanup():
    # Before test section
    yield
    # After test section
    _cleanup()


@pytest.fixture
def tmp_dir(tmpdir):
    # tmpdir is py.path.LocalPath
    # * python 3.6 supports non string/bytes object as path in os.path.*
    #   methods
    #
    # * python 3.5 fails if gets non string/bytes object in os.path.* methods
    #   E.g. you will get error like:
    #     TypeError: join() argument must be str or bytes, not 'LocalPath'
    #
    # So in order to support 3.5, let's create fixture which returns string obj
    return str(tmpdir)


@pytest.fixture
def tmp_path(tmpdir):
    return Path(str(tmpdir))
