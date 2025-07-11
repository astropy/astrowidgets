import pytest

from traitlets import TraitError

from astro_image_display_api import ImageAPITest
from astro_image_display_api import ImageViewerInterface

_ = pytest.importorskip("bqplot",
                        reason="Package required for test is not "
                               "available.")
from astrowidgets.bqplot import ImageWidget  # noqa: E402

def test_instance():
    image = ImageWidget()
    assert isinstance(image, ImageViewerInterface)

class TestBQplotWidget(ImageAPITest):
    image_widget_class = ImageWidget
    cursor_error_classes = (ValueError, TraitError)

    @pytest.mark.skip(reason="Saving is done in javascript and requires "
                             "a running browser.")
    def test_save(self, tmp_path):
        pass

    @pytest.mark.skip(reason="Saving is done in javascript and requires "
                             "a running browser.")
    def test_save_overwrite(self, tmp_path):
        pass
