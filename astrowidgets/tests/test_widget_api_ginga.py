import pytest

from .widget_api_test import ImageWidgetAPITest
from astrowidgets.interface_definition import ImageViewerInterface

ginga = pytest.importorskip("ginga",
                            reason="Package required for test is not "
                                   "available.")
from astrowidgets.ginga import ImageWidget  # noqa: E402

def test_instance():
    image = ImageWidget()
    assert isinstance(image, ImageViewerInterface)

class TestGingaWidget(ImageWidgetAPITest):
    image_widget_class = ImageWidget
