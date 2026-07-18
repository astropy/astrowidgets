import astropy.visualization as apviz
from astropy.table import Table
import numpy as np
import pytest

from traitlets import TraitError

from astro_image_display_api.api_test import ImageAPITest
from astro_image_display_api import ImageViewerInterface

bqplot = pytest.importorskip("bqplot",
                             reason="Package required for test is not "
                                    "available.")
from astrowidgets.bqplot import ImageWidget, bqcolors  # noqa: E402


def test_instance():
    image = ImageWidget()
    assert isinstance(image, ImageViewerInterface)


def test_mouse_click_does_not_raise_or_block_callbacks():
    # Regression test for #206: the built-in click handler referenced
    # attributes that no longer exist, so any click raised AttributeError
    # and, because ipywidgets runs on_msg callbacks in registration order
    # with no exception isolation, blocked user-registered callbacks.
    image = ImageWidget()

    # A click before any image is loaded should be a no-op.
    image._mouse_click({'domain': {'x': 3, 'y': 3}})

    image.load_image(np.zeros((10, 10)))

    calls = []

    def user_callback(interaction, event_data, buffers):
        calls.append(event_data)

    image._astro_im.interaction.on_msg(user_callback)

    # Simulate the comm message a real mouse click produces; this invokes
    # all registered on_msg callbacks in order, built-in handler first.
    click_event = {'event': 'click', 'domain': {'x': 3, 'y': 3}}
    image._astro_im.interaction._handle_custom_msg(click_event, [])

    assert calls == [click_event]


@pytest.mark.parametrize(
    "shape", ["circle", "square", "crosshair", "plus", "diamond"]
)
def test_catalog_markers_use_scatter_shape_size_and_arrays(shape):
    image = ImageWidget()

    image.load_catalog(
        Table({"x": [1.0], "y": [2.0]}),
        catalog_label="test",
        catalog_style={"color": "red", "shape": shape, "size": 6},
    )
    marker = image._astro_im._scatter_marks["test"]
    assert type(marker) is bqplot.Scatter
    assert marker.default_size == 36
    np.testing.assert_array_equal(marker.x, [1.0])
    np.testing.assert_array_equal(marker.y, [2.0])

    image.set_catalog_style(
        catalog_label="test", size=7, shape=shape
    )

    marker = image._astro_im._scatter_marks["test"]
    assert type(marker) is bqplot.Scatter
    assert marker.default_size == 49
    assert marker.marker == shape


def test_catalog_marks_use_resolved_label():
    # Regression test for #213: set_catalog_style and remove_catalog used
    # the caller's (unresolved) catalog_label as the bqplot mark id, so an
    # unlabeled set_catalog_style call plotted an orphan second scatter
    # under the id "None" and an unlabeled remove_catalog raised. The
    # resolved label must always be the mark id.
    image = ImageWidget()
    image.load_catalog(Table({'x': [1.0], 'y': [2.0]}), catalog_label='test')
    assert list(image._astro_im._scatter_marks) == ['test']

    image.set_catalog_style(size=7)
    assert list(image._astro_im._scatter_marks) == ['test']
    assert image._astro_im._scatter_marks['test'].default_size == 49

    image.remove_catalog()
    assert image._astro_im._scatter_marks == {}


def test_catalog_generated_label_targets_marks():
    # A catalog loaded without a label gets a generated label; styling or
    # removing it by that generated label must target its marks.
    image = ImageWidget()
    image.load_catalog(Table({'x': [1.0], 'y': [2.0]}))

    label = image.catalog_labels[0]
    assert list(image._astro_im._scatter_marks) == [label]

    image.set_catalog_style(catalog_label=label, size=7)
    assert list(image._astro_im._scatter_marks) == [label]
    assert image._astro_im._scatter_marks[label].default_size == 49

    image.remove_catalog(catalog_label=label)
    assert image._astro_im._scatter_marks == {}


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
        # both before and after an image is loaded. Settings are stored per
        # image, so get_colormap can only report it once an image is loaded.
        greys_r = bqcolors('Greys_r')
        assert self.image._astro_im._image.scales['image'].colors == greys_r

        self.image.load_image(data)
        assert self.image.get_colormap() == 'Greys_r'
        assert self.image._astro_im._image.scales['image'].colors == greys_r

    def test_load_image_new_label_gets_default_settings(self):
        # Loading an image under a new label should display it with the
        # widget's default cuts, stretch and colormap -- not settings
        # carried forward from the previously displayed image -- and the
        # previous label keeps its own settings untouched.
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

        second_cuts = self.image.get_cuts(image_label='second')
        assert isinstance(second_cuts, apviz.AsymmetricPercentileInterval)
        assert second_cuts.lower_percentile == 30
        assert second_cuts.upper_percentile == 96
        assert isinstance(self.image.get_stretch(image_label='second'),
                          apviz.LinearStretch)
        assert self.image.get_colormap(image_label='second') == 'Greys_r'
        assert self.image._astro_im._image.scales['image'].colors == bqcolors('Greys_r')

        displayed = np.asarray(self.image._astro_im._image.image)
        np.testing.assert_allclose(displayed, second_cuts(arr))

        # The first label's settings are untouched.
        assert self.image.get_cuts(image_label='first') is cuts
        assert self.image.get_stretch(image_label='first') is stretch
        assert self.image.get_colormap(image_label='first') == 'viridis'

    def test_reload_existing_label_keeps_its_settings(self):
        # Loading new data under an existing image label keeps the settings
        # already stored for that label -- not the settings of the image it
        # replaces on the display -- and redisplays with them.
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
        # every image gets a generated label. Loading a second image must
        # display it with the settings currently in effect, not fall back
        # to the widget defaults.
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

        # Both images stay loaded under generated labels, so the get_*
        # methods need a label; ask about the displayed (new) image.
        displayed_label = self.image._displayed_image_labels[0]
        assert self.image.get_cuts(image_label=displayed_label) is cuts
        assert self.image.get_stretch(image_label=displayed_label) is stretch
        assert self.image.get_colormap(image_label=displayed_label) == 'viridis'
        assert self.image._astro_im._image.scales['image'].colors == bqcolors('viridis')

        displayed = np.asarray(self.image._astro_im._image.image)
        np.testing.assert_allclose(displayed, stretch(cuts(arr2)))

    def test_load_image_without_label_never_replaces(self):
        # An unlabeled load gets a generated label, so it can never
        # silently replace a previously loaded image: the labeled image
        # keeps its data and settings, while the new image takes over the
        # display with the widget's default settings.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        labeled_cuts = apviz.ManualInterval(1100, 1300)
        self.image.load_image(arr, image_label='labeled')
        self.image.set_cuts(labeled_cuts, image_label='labeled')

        self.image.load_image(arr + 2)

        # The labeled image is untouched...
        np.testing.assert_array_equal(
            np.asarray(self.image.get_image(image_label='labeled')), arr)
        assert self.image.get_cuts(image_label='labeled') is labeled_cuts

        # ...and the new image is displayed, under a generated label, with
        # the widget's default cuts rather than the labeled image's.
        new_label = self.image._displayed_image_labels[0]
        assert new_label != 'labeled'
        np.testing.assert_array_equal(
            np.asarray(self.image.get_image(image_label=new_label)), arr + 2)
        new_cuts = self.image.get_cuts(image_label=new_label)
        assert new_cuts is not labeled_cuts
        assert isinstance(new_cuts, apviz.AsymmetricPercentileInterval)
        assert new_cuts.lower_percentile == 30
        assert new_cuts.upper_percentile == 96

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
        # Cuts that differ from the widget default, so that a refresh
        # falling back to the default cuts produces a different array.
        cuts = apviz.AsymmetricPercentileInterval(5, 90)
        self.image.set_cuts(cuts)
        stretch = apviz.LogStretch()
        self.image.set_stretch(stretch)

        displayed = np.asarray(self.image._astro_im._image.image)
        np.testing.assert_allclose(displayed, stretch(cuts(arr)))

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

        # Batching must not change the end state: the new label is displayed
        # with the widget's default cuts.
        cuts = self.image.get_cuts(image_label='second')
        assert isinstance(cuts, apviz.AsymmetricPercentileInterval)
        displayed = np.asarray(image_mark.image)
        np.testing.assert_allclose(displayed, cuts(arr))
        np.testing.assert_allclose(image_mark.x, [-0.5, arr.shape[1] - 0.5])
        np.testing.assert_allclose(image_mark.y, [-0.5, arr.shape[0] - 0.5])

    def test_load_image_flushes_image_before_scales(self, data, mocker):
        # The image mark and the two scales are separate widgets, each
        # syncing its own state message that the front end redraws on.
        # If a scale syncs before the image mark, the front end briefly
        # draws the OLD image against the NEW scales (a visible refit).
        # The new image data must reach the front end first.
        self.image.load_image(data, image_label='first')

        astro_im = self.image._astro_im
        image_mark = astro_im._image
        scale_x = astro_im._scales['x']
        scale_y = astro_im._scales['y']

        # A different shape so the image extent and both scales change.
        arr = np.arange(30 * 40, dtype=float).reshape(30, 40)

        order = []

        def record(name, widget):
            real = widget.send_state

            def wrapper(*args, **kwargs):
                order.append(name)
                return real(*args, **kwargs)

            mocker.patch.object(widget, 'send_state', side_effect=wrapper)

        record('image', image_mark)
        record('scale_x', scale_x)
        record('scale_y', scale_y)

        self.image.load_image(arr, image_label='second')

        assert 'image' in order
        assert order.index('image') < order.index('scale_x')
        assert order.index('image') < order.index('scale_y')

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
