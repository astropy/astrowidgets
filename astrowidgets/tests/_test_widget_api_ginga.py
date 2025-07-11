import pytest

from astro_image_display_api import ImageWidgetAPITest
from astro_image_display_api import ImageViewerInterface

ginga = pytest.importorskip("ginga",
                            reason="Package required for test is not "
                                   "available.")
from astrowidgets.ginga import ImageWidget  # noqa: E402


def test_instance():
    image = ImageWidget()
    assert isinstance(image, ImageViewerInterface)


class TestGingaWidget(ImageWidgetAPITest):
    image_widget_class = ImageWidget
