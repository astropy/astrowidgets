import numpy as np
import pytest

from traitlets import TraitError

from astro_image_display_api.api_test import ImageAPITest
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

    def test_initial_display_uses_stored_cuts(self):
        # A high-dynamic-range image: narrow uint16 background plus a single
        # saturated pixel. Scaling to the full data range would render the
        # background essentially black, so the initial display must use the
        # stored cuts (a percentile interval), not min/max.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        self.image.load_image(arr)

        displayed = np.asarray(self.image._astro_im._image.image)
        expected = self.image.get_cuts()(arr)
        np.testing.assert_allclose(displayed, expected)

    def test_get_viewport_reflects_interactive_pan(self, data):
        # Panning in the browser shifts the bqplot scales directly. Simulate
        # that here and check that get_viewport reports the new center.
        self.image.load_image(data)
        self.image.set_viewport(center=(75, 50), fov=100)

        scales = self.image._astro_im._scales

        # Horizontal pan: shift the x scale. Set min before max so that the
        # center is already correct when the 'max' observer fires.
        dx = 20
        scales['x'].min += dx
        scales['x'].max += dx

        vp = self.image.get_viewport(sky_or_pixel='pixel')
        assert vp['center'][0] == pytest.approx(75 + dx)
        assert vp['center'][1] == pytest.approx(50)

        # Vertical pan: shift the y scale.
        dy = -15
        scales['y'].min += dy
        scales['y'].max += dy

        vp = self.image.get_viewport(sky_or_pixel='pixel')
        assert vp['center'][0] == pytest.approx(75 + dx)
        assert vp['center'][1] == pytest.approx(50 + dy)
