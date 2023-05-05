import pytest

from traitlets import TraitError

from .widget_api_test import ImageWidgetAPITest

_ = pytest.importorskip("bqplot",
                        reason="Package required for test is not "
                               "available.")
from astrowidgets.bqplot import ImageWidget  # noqa: E402


class TestBQplotWidget(ImageWidgetAPITest):
    image_widget_class = ImageWidget
    cursor_error_classes = (ValueError, TraitError)

    @pytest.mark.skip(reason="Saving is done in javascript and requires "
                             "a running browser.")
    def test_save(self, tmpdir):
        pass
