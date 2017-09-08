import os
import pytest

from bed.bedtrace import cleanup


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
    cleanup()
