import pytest
import downstream.aging as aging


def test_is_od():
    assert aging.is_od("FOO:OD1") is True


def test_is_yd():
    assert aging.is_yd("FOO:YD1") is True


@pytest.mark.parametrize("value,expected", [
    ("foo.od17.boo", "od17"),
    ("foo.OD17.boo", "OD17"),
    ("foo.YD17.boo", "YD17"),
    ("foo.OD_OD5.boo", "OD5"),
    ("foo.ODS.boo", None),
])
def test_age(value, expected):
    assert aging.donor(value) == expected
