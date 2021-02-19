"""Test to make sure astrowidgets can install and be used with
only abstract class, without optional backend.

"""
import pytest

from astrowidgets.core import BaseImageWidget


class DummyWidget(BaseImageWidget):
    pass


def test_abstract_no_imp():
    # Ensure subclass cannot be used without implementing all the abstracted
    # things.
    with pytest.raises(TypeError, match="Can't instantiate abstract class"):
        DummyWidget()
