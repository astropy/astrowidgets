import astropy.visualization as apviz
import numpy as np
import pytest

from traitlets import TraitError

from astro_image_display_api.api_test import ImageAPITest
from astro_image_display_api import ImageViewerInterface

_ = pytest.importorskip("bqplot",
                        reason="Package required for test is not "
                               "available.")
from astrowidgets.bqplot import ImageWidget, bqcolors  # noqa: E402


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

    def test_default_cuts_are_30_96_percentile(self):
        # When no cuts are stored or passed, the widget's fallback should
        # be the 30-96 percentile interval, which cuts out the sky
        # background and clips only the brightest pixels.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        self.image._data = arr
        expected = apviz.AsymmetricPercentileInterval(30, 96)(arr)
        np.testing.assert_allclose(self.image._interval_and_stretch(), expected)

    def test_default_colormap_is_greys_r(self, data):
        # With no colormap explicitly set, the display should use Greys_r
        # (low = black, high = white), the usual astronomical convention,
        # both before and after an image is loaded, and get_colormap should
        # report it.
        greys_r = bqcolors('Greys_r')
        assert self.image._astro_im._image.scales['image'].colors == greys_r

        self.image.load_image(data)
        assert self.image.get_colormap() == 'Greys_r'
        assert self.image._astro_im._image.scales['image'].colors == greys_r

    def test_load_image_keeps_current_display_settings(self):
        # Loading a new image should display it with the cuts, stretch and
        # colormap that are currently in effect, carried forward from the
        # previously displayed image, and store them for the new image so
        # that the get_* methods agree with the display. Only the viewport
        # resets on load.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        self.image.load_image(arr, image_label='first')
        cuts = apviz.AsymmetricPercentileInterval(5, 90)
        stretch = apviz.LogStretch()
        self.image.set_cuts(cuts, image_label='first')
        self.image.set_stretch(stretch, image_label='first')
        self.image.set_colormap('viridis', image_label='first')

        self.image.load_image(arr, image_label='second')

        assert self.image.get_cuts(image_label='second') is cuts
        assert self.image.get_stretch(image_label='second') is stretch
        assert self.image.get_colormap(image_label='second') == 'viridis'
        assert self.image._astro_im._image.scales['image'].colors == bqcolors('viridis')

        displayed = np.asarray(self.image._astro_im._image.image)
        np.testing.assert_allclose(displayed, stretch(cuts(arr)))

    def test_reload_existing_label_keeps_its_settings(self):
        # When new data is loaded under an existing image label, that
        # label's own stored settings are the current ones for the image
        # and must be kept, even if a different image was loaded (and so
        # displayed) more recently.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        first_cuts = apviz.AsymmetricPercentileInterval(5, 90)
        first_stretch = apviz.LogStretch()
        self.image.load_image(arr, image_label='first')
        self.image.set_cuts(first_cuts, image_label='first')
        self.image.set_stretch(first_stretch, image_label='first')
        self.image.set_colormap('viridis', image_label='first')

        self.image.load_image(arr, image_label='second')
        self.image.set_cuts(apviz.ManualInterval(1100, 1300),
                            image_label='second')
        self.image.set_stretch(apviz.SqrtStretch(), image_label='second')
        self.image.set_colormap('plasma', image_label='second')

        new_arr = arr + 10
        self.image.load_image(new_arr, image_label='first')

        assert self.image.get_cuts(image_label='first') is first_cuts
        assert self.image.get_stretch(image_label='first') is first_stretch
        assert self.image.get_colormap(image_label='first') == 'viridis'
        assert self.image._astro_im._image.scales['image'].colors == bqcolors('viridis')

        displayed = np.asarray(self.image._astro_im._image.image)
        np.testing.assert_allclose(displayed, first_stretch(first_cuts(new_arr)))

    def test_load_image_keeps_settings_without_labels(self):
        # The carry-forward of cuts, stretch and colormap must work in the
        # common interactive case where no image_label is ever passed, so
        # every image and its settings live under the API's default (None)
        # label. Loading a second image must keep the settings currently in
        # effect, not fall back to the widget defaults.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        self.image.load_image(arr)
        cuts = apviz.AsymmetricPercentileInterval(5, 90)
        stretch = apviz.LogStretch()
        self.image.set_cuts(cuts)
        self.image.set_stretch(stretch)
        self.image.set_colormap('viridis')

        # A differently shaped second image, still with no label.
        arr2 = rng.integers(1100, 1300, size=(40, 30)).astype(np.uint16)
        arr2[10, 15] = 65535
        self.image.load_image(arr2)

        assert self.image.get_cuts() is cuts
        assert self.image.get_stretch() is stretch
        assert self.image.get_colormap() == 'viridis'
        assert self.image._astro_im._image.scales['image'].colors == bqcolors('viridis')

        displayed = np.asarray(self.image._astro_im._image.image)
        np.testing.assert_allclose(displayed, stretch(cuts(arr2)))

    def test_first_load_stores_widget_default_cuts(self):
        # With nothing loaded yet there are no current settings to carry
        # forward, so the first image should be displayed with the widget's
        # default cuts, and those cuts should be stored so that get_cuts
        # reports what is displayed.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        self.image.load_image(arr)

        cuts = self.image.get_cuts()
        assert isinstance(cuts, apviz.AsymmetricPercentileInterval)
        assert cuts.lower_percentile == 30
        assert cuts.upper_percentile == 96

        displayed = np.asarray(self.image._astro_im._image.image)
        np.testing.assert_allclose(displayed, cuts(arr))

    def test_set_stretch_keeps_stored_cuts(self):
        # Changing the stretch must re-display using the cuts stored for
        # the image, not fall back to the widget's default cuts.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        self.image.load_image(arr)
        stretch = apviz.LogStretch()
        self.image.set_stretch(stretch)

        displayed = np.asarray(self.image._astro_im._image.image)
        expected = stretch(self.image.get_cuts()(arr))
        np.testing.assert_allclose(displayed, expected)

    def test_set_cuts_keeps_stored_stretch(self):
        # Changing the cuts must re-display using the stretch stored for
        # the image, not fall back to the widget's default stretch.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        self.image.load_image(arr)
        stretch = apviz.LogStretch()
        self.image.set_stretch(stretch)
        cuts = apviz.AsymmetricPercentileInterval(5, 90)
        self.image.set_cuts(cuts)

        displayed = np.asarray(self.image._astro_im._image.image)
        expected = stretch(cuts(arr))
        np.testing.assert_allclose(displayed, expected)

    def test_load_image_batches_frontend_updates(self, data, mocker):
        # Every traitlet assignment is normally synced to the browser as
        # its own message, each triggering a redraw, so loading an image
        # produced a series of visible intermediate states (flicker).
        # Loading should batch the updates so the image mark and each
        # scale send at most one state message.
        self.image.load_image(data, image_label='first')

        astro_im = self.image._astro_im
        image_mark = astro_im._image
        scale_x = astro_im._scales['x']
        scale_y = astro_im._scales['y']

        # Use a different shape so the image extent and scales all change.
        arr = np.arange(30 * 40, dtype=float).reshape(30, 40)

        spy_image = mocker.spy(image_mark, 'send_state')
        spy_x = mocker.spy(scale_x, 'send_state')
        spy_y = mocker.spy(scale_y, 'send_state')

        self.image.load_image(arr, image_label='second')

        assert spy_image.call_count <= 1
        assert spy_x.call_count <= 1
        assert spy_y.call_count <= 1

        # Batching must not change the end state.
        displayed = np.asarray(image_mark.image)
        expected = self.image.get_cuts(image_label='second')(arr)
        np.testing.assert_allclose(displayed, expected)
        np.testing.assert_allclose(image_mark.x, [-0.5, arr.shape[1] - 0.5])
        np.testing.assert_allclose(image_mark.y, [-0.5, arr.shape[0] - 0.5])

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
