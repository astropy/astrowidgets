import pytest

import ipywidgets as ipyw
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


def test_widget_kwargs_forwarded_to_vbox():
    # Widget keyword arguments should be passed through to the underlying
    # ipywidgets.VBox, so e.g. a custom layout is honored.
    image = ImageWidget(layout=ipyw.Layout(border='1px solid red'))
    assert image.layout.border == '1px solid red'


def test_use_opencv_kwarg_warns_not_errors():
    # The deprecated use_opencv kwarg should warn rather than be forwarded to
    # VBox (which would raise a TraitError for an unknown trait).
    with pytest.warns(DeprecationWarning):
        ImageWidget(use_opencv=True)


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
