import pytest

from traitlets import TraitError

from astro_image_display_api.api_test import ImageAPITest
from astro_image_display_api import ImageViewerInterface

_ = pytest.importorskip("ginga",
                        reason="Package required for test is not "
                               "available.")
from astrowidgets.ginga import ImageWidget  # noqa: E402


def test_instance():
    image = ImageWidget()
    assert isinstance(image, ImageViewerInterface)


class TestGingaWidget(ImageAPITest):
    image_widget_class = ImageWidget
    cursor_error_classes = (ValueError, TraitError)

    @pytest.mark.skip(reason="Saving requires a running browser to populate "
                             "the image buffer.")
    def test_save(self, tmp_path):
        pass

    @pytest.mark.skip(reason="Saving requires a running browser to populate "
                             "the image buffer.")
    def test_save_overwrite(self, tmp_path):
        pass
