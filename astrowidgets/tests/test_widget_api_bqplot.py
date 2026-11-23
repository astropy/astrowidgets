import astropy.visualization as apviz
import numpy as np
import pytest
from astro_image_display_api import ImageViewerInterface
from astro_image_display_api.api_test import ImageAPITest
from astropy.table import Table
from matplotlib.colors import to_rgb
from traitlets import TraitError

bqplot = pytest.importorskip(
    "bqplot", reason="Package required for test is not available."
)
from astrowidgets.bqplot import (
    ImageWidget,
    bqcolors,
    scale_to_levels,
    stretch_fits_in_colors,
)


def assert_displayed(widget, arr, cuts, stretch=None, colormap="Greys_r"):
    """
    Check what a widget has sent the front end to draw an image with.

    The front end draws a pixel in color number
    ``(value - min) / (max - min)`` of the color scale, so what is on the
    screen is fixed by the pixels, the limits and the colors it is sent.
    Which of the two ways the widget sends them depends on the stretch.

    Parameters
    ----------
    widget : `astrowidgets.bqplot.ImageWidget`
        The widget to check.
    arr : array-like
        The image that should be displayed.
    cuts : `astropy.visualization.BaseInterval`
        The cuts whose limits the color scale should run between.
    stretch : `astropy.visualization.BaseStretch`, optional
        The stretch that should be built into the colors.
    colormap : str, optional
        The colormap the colors should be from.
    """
    image_mark = widget._astro_im._image
    color_scale = image_mark.scales["image"]
    arr = np.asarray(arr)
    low, high = cuts.get_limits(arr)

    if stretch is None or stretch_fits_in_colors(stretch):
        # The pixels as they are, the cuts as the limits of the color
        # scale and the stretch as the spacing of its colors.
        assert image_mark.image.dtype == np.float32
        np.testing.assert_array_equal(image_mark.image, arr)
        np.testing.assert_allclose(
            (color_scale.min, color_scale.max), (low, high)
        )
        assert color_scale.colors == bqcolors(colormap, stretch=stretch)
    else:
        # The pixels with the cuts and stretch applied, as the number of
        # the color each is drawn in, and evenly spaced colors.
        assert image_mark.image.dtype == np.uint8
        scaled = stretch(apviz.ManualInterval(low, high)(arr))
        np.testing.assert_array_equal(
            image_mark.image, np.minimum(np.floor(scaled * 256), 255)
        )
        assert (color_scale.min, color_scale.max) == (0, 1)
        assert color_scale.colors == bqcolors(colormap)


def screen_levels(widget):
    """
    Work out what the front end draws, as the fragment shader does.

    Parameters
    ----------
    widget : `astrowidgets.bqplot.ImageWidget`
        A widget showing an image in the ``gray`` colormap, whose
        brightness rises evenly from black to white. (``Greys_r``, the
        default, does not rise evenly.)

    Returns
    -------
    `numpy.ndarray`
        The brightness each pixel is drawn at, from 0 to 1.
    """
    image_mark = widget._astro_im._image
    color_scale = image_mark.scales["image"]
    pixels = np.asarray(image_mark.image)
    if pixels.dtype == np.uint8:
        # The GPU reads an 8 bit texture as a fraction of 255.
        pixels = pixels / 255
    value = (pixels - color_scale.min) / (color_scale.max - color_scale.min)
    # The colors are a texture 256 wide, read at the nearest texel and
    # clamped at its ends.
    index = np.clip(np.floor(value * 256), 0, 255).astype(int)
    greys = np.array([to_rgb(color)[0] for color in color_scale.colors])
    return greys[index]


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
    image._mouse_click({"domain": {"x": 3, "y": 3}})

    image.load_image(np.zeros((10, 10)))

    calls = []

    def user_callback(interaction, event_data, buffers):
        calls.append(event_data)

    image._astro_im.interaction.on_msg(user_callback)

    # Simulate the comm message a real mouse click produces; this invokes
    # all registered on_msg callbacks in order, built-in handler first.
    click_event = {"event": "click", "domain": {"x": 3, "y": 3}}
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

    image.set_catalog_style(catalog_label="test", size=7, shape=shape)

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
    image.load_catalog(Table({"x": [1.0], "y": [2.0]}), catalog_label="test")
    assert list(image._astro_im._scatter_marks) == ["test"]

    image.set_catalog_style(size=7)
    assert list(image._astro_im._scatter_marks) == ["test"]
    assert image._astro_im._scatter_marks["test"].default_size == 49

    image.remove_catalog()
    assert image._astro_im._scatter_marks == {}


def test_catalog_default_label_targets_marks():
    # A catalog loaded without a label goes under the shared default
    # label; styling or removing it by that label must target its marks.
    image = ImageWidget()
    image.load_catalog(Table({"x": [1.0], "y": [2.0]}))

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

    @pytest.mark.skip(
        reason="Saving is done in javascript and requires a running browser."
    )
    def test_save(self, tmp_path):
        pass

    @pytest.mark.skip(
        reason="Saving is done in javascript and requires a running browser."
    )
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

        assert_displayed(self.image, arr, self.image.get_cuts())
        # The background is spread over the colors, not crushed into the
        # first few of them by the saturated pixel.
        self.image.set_colormap("gray")
        assert np.median(screen_levels(self.image)) > 0.2

    def test_default_cuts_are_30_96_percentile(self):
        # When no cuts are stored or passed, the widget's fallback should
        # be the 30-96 percentile interval, which cuts out the sky
        # background and clips only the brightest pixels.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        self.image.load_image(arr)
        assert_displayed(
            self.image, arr, apviz.AsymmetricPercentileInterval(30, 96)
        )

    def test_default_colormap_is_greys_r(self, data):
        # With no colormap explicitly set, the display should use Greys_r
        # (low = black, high = white), the usual astronomical convention,
        # both before and after an image is loaded. Settings are stored per
        # image, so get_colormap can only report it once an image is loaded.
        greys_r = bqcolors("Greys_r")
        assert self.image._astro_im._image.scales["image"].colors == greys_r

        self.image.load_image(data)
        assert self.image.get_colormap() == "Greys_r"
        assert self.image._astro_im._image.scales["image"].colors == greys_r

    def test_load_image_new_label_gets_default_settings(self):
        # Loading an image under a new label should display it with the
        # widget's default cuts, stretch and colormap -- not settings
        # carried forward from the previously displayed image -- and the
        # previous label keeps its own settings untouched.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        self.image.load_image(arr, image_label="first")
        cuts = apviz.AsymmetricPercentileInterval(5, 90)
        stretch = apviz.LogStretch()
        self.image.set_cuts(cuts, image_label="first")
        self.image.set_stretch(stretch, image_label="first")
        self.image.set_colormap("viridis", image_label="first")

        self.image.load_image(arr, image_label="second")

        second_cuts = self.image.get_cuts(image_label="second")
        assert isinstance(second_cuts, apviz.AsymmetricPercentileInterval)
        assert second_cuts.lower_percentile == 30
        assert second_cuts.upper_percentile == 96
        assert isinstance(
            self.image.get_stretch(image_label="second"), apviz.LinearStretch
        )
        assert self.image.get_colormap(image_label="second") == "Greys_r"
        assert self.image._astro_im._image.scales["image"].colors == bqcolors(
            "Greys_r"
        )

        assert_displayed(
            self.image, arr, second_cuts, stretch=apviz.LinearStretch()
        )

        # The first label's settings are untouched.
        assert self.image.get_cuts(image_label="first") is cuts
        assert self.image.get_stretch(image_label="first") is stretch
        assert self.image.get_colormap(image_label="first") == "viridis"

    def test_reload_existing_label_keeps_its_settings(self):
        # Loading new data under an existing image label keeps the settings
        # already stored for that label -- not the settings of the image it
        # replaces on the display -- and redisplays with them.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        first_cuts = apviz.AsymmetricPercentileInterval(5, 90)
        first_stretch = apviz.LogStretch()
        self.image.load_image(arr, image_label="first")
        self.image.set_cuts(first_cuts, image_label="first")
        self.image.set_stretch(first_stretch, image_label="first")
        self.image.set_colormap("viridis", image_label="first")

        self.image.load_image(arr, image_label="second")
        self.image.set_cuts(
            apviz.ManualInterval(1100, 1300), image_label="second"
        )
        self.image.set_stretch(apviz.SqrtStretch(), image_label="second")
        self.image.set_colormap("plasma", image_label="second")

        new_arr = arr + 10
        self.image.load_image(new_arr, image_label="first")

        assert self.image.get_cuts(image_label="first") is first_cuts
        assert self.image.get_stretch(image_label="first") is first_stretch
        assert self.image.get_colormap(image_label="first") == "viridis"
        assert_displayed(
            self.image,
            new_arr,
            first_cuts,
            stretch=first_stretch,
            colormap="viridis",
        )

    def test_load_image_keeps_settings_without_labels(self):
        # The carry-forward of cuts, stretch and colormap must work in the
        # common interactive case where no image_label is ever passed, so
        # every load targets the shared default label. Loading a second
        # image must display it with the settings currently in effect, not
        # fall back to the widget defaults.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        self.image.load_image(arr)
        cuts = apviz.AsymmetricPercentileInterval(5, 90)
        stretch = apviz.LogStretch()
        self.image.set_cuts(cuts)
        self.image.set_stretch(stretch)
        self.image.set_colormap("viridis")

        # A differently shaped second image, still with no label.
        arr2 = rng.integers(1100, 1300, size=(40, 30)).astype(np.uint16)
        arr2[10, 15] = 65535
        self.image.load_image(arr2)

        # The second load replaced the first under the shared default
        # label, keeping the settings attached to that label.
        displayed_label = self.image._displayed_image_labels[0]
        assert self.image.get_cuts(image_label=displayed_label) is cuts
        assert self.image.get_stretch(image_label=displayed_label) is stretch
        assert (
            self.image.get_colormap(image_label=displayed_label) == "viridis"
        )
        assert_displayed(
            self.image, arr2, cuts, stretch=stretch, colormap="viridis"
        )

    def test_load_image_without_label_targets_default_label(self):
        # An unlabeled load resolves to the shared default label, so it
        # never silently replaces an explicitly labeled image: the labeled
        # image keeps its data and settings, while the new image takes
        # over the display with the widget's default settings. Repeated
        # unlabeled loads, though, replace each other under that one
        # default label.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        labeled_cuts = apviz.ManualInterval(1100, 1300)
        self.image.load_image(arr, image_label="labeled")
        self.image.set_cuts(labeled_cuts, image_label="labeled")

        self.image.load_image(arr + 2)

        # The labeled image is untouched...
        np.testing.assert_array_equal(
            np.asarray(self.image.get_image(image_label="labeled")), arr
        )
        assert self.image.get_cuts(image_label="labeled") is labeled_cuts

        # ...and the new image is displayed, under the default label, with
        # the widget's default cuts rather than the labeled image's.
        new_label = self.image._displayed_image_labels[0]
        assert new_label != "labeled"
        np.testing.assert_array_equal(
            np.asarray(self.image.get_image(image_label=new_label)), arr + 2
        )
        new_cuts = self.image.get_cuts(image_label=new_label)
        assert new_cuts is not labeled_cuts
        assert isinstance(new_cuts, apviz.AsymmetricPercentileInterval)
        assert new_cuts.lower_percentile == 30
        assert new_cuts.upper_percentile == 96

        # A further unlabeled load replaces the previous unlabeled image
        # in place: same label, new data, no third image.
        self.image.load_image(arr + 5)
        assert self.image._displayed_image_labels[0] == new_label
        np.testing.assert_array_equal(
            np.asarray(self.image.get_image(image_label=new_label)), arr + 5
        )
        assert sorted(self.image.image_labels) == sorted(
            ("labeled", new_label)
        )

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

        assert_displayed(self.image, arr, cuts)

    def test_set_stretch_keeps_stored_cuts(self):
        # Changing the stretch must leave the display using the cuts stored
        # for the image, not fall back to the widget's default cuts.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        self.image.load_image(arr)
        # Cuts that differ from the widget default, so that a refresh
        # falling back to the default cuts gives different color limits.
        cuts = apviz.AsymmetricPercentileInterval(5, 90)
        self.image.set_cuts(cuts)
        stretch = apviz.LogStretch()
        self.image.set_stretch(stretch)

        assert_displayed(self.image, arr, cuts, stretch=stretch)

    def test_set_cuts_keeps_stored_stretch(self):
        # Changing the cuts must leave the display using the stretch stored
        # for the image, not fall back to the widget's default stretch.
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        arr[25, 30] = 65535

        self.image.load_image(arr)
        stretch = apviz.LogStretch()
        self.image.set_stretch(stretch)
        cuts = apviz.AsymmetricPercentileInterval(5, 90)
        self.image.set_cuts(cuts)

        assert_displayed(self.image, arr, cuts, stretch=stretch)

    def test_load_image_batches_frontend_updates(self, data, mocker):
        # Every traitlet assignment is normally synced to the browser as
        # its own message, each triggering a redraw, so loading an image
        # produced a series of visible intermediate states (flicker).
        # Loading should batch the updates so the image mark and each
        # scale send at most one state message.
        self.image.load_image(data, image_label="first")

        astro_im = self.image._astro_im
        image_mark = astro_im._image
        scale_x = astro_im._scales["x"]
        scale_y = astro_im._scales["y"]

        # Use a different shape so the image extent and scales all change.
        arr = np.arange(30 * 40, dtype=float).reshape(30, 40)

        spy_image = mocker.spy(image_mark, "send_state")
        spy_x = mocker.spy(scale_x, "send_state")
        spy_y = mocker.spy(scale_y, "send_state")

        self.image.load_image(arr, image_label="second")

        assert spy_image.call_count <= 1
        assert spy_x.call_count <= 1
        assert spy_y.call_count <= 1

        # Batching must not change the end state: the new label is displayed
        # with the widget's default cuts.
        cuts = self.image.get_cuts(image_label="second")
        assert isinstance(cuts, apviz.AsymmetricPercentileInterval)
        assert_displayed(self.image, arr, cuts)
        np.testing.assert_allclose(image_mark.x, [-0.5, arr.shape[1] - 0.5])
        np.testing.assert_allclose(image_mark.y, [-0.5, arr.shape[0] - 0.5])

    def test_load_image_flushes_image_before_scales(self, data, mocker):
        # The image mark and the two scales are separate widgets, each
        # syncing its own state message that the front end redraws on.
        # If a scale syncs before the image mark, the front end briefly
        # draws the OLD image against the NEW scales (a visible refit).
        # The new image data must reach the front end first.
        self.image.load_image(data, image_label="first")

        astro_im = self.image._astro_im
        image_mark = astro_im._image
        scale_x = astro_im._scales["x"]
        scale_y = astro_im._scales["y"]

        # A different shape so the image extent and both scales change.
        arr = np.arange(30 * 40, dtype=float).reshape(30, 40)

        order = []

        def record(name, widget):
            real = widget.send_state

            def wrapper(*args, **kwargs):
                order.append(name)
                return real(*args, **kwargs)

            mocker.patch.object(widget, "send_state", side_effect=wrapper)

        record("image", image_mark)
        record("scale_x", scale_x)
        record("scale_y", scale_y)

        self.image.load_image(arr, image_label="second")

        assert "image" in order
        assert order.index("image") < order.index("scale_x")
        assert order.index("image") < order.index("scale_y")

    @pytest.mark.parametrize(
        "change",
        [
            lambda image: image.set_cuts((1150, 1250)),
            lambda image: image.set_cuts(
                apviz.AsymmetricPercentileInterval(5, 90)
            ),
            lambda image: image.set_stretch(apviz.SqrtStretch()),
            lambda image: image.set_colormap("viridis"),
            lambda image: image.remove_catalog(catalog_label="*"),
        ],
        ids=[
            "manual cuts",
            "percentile cuts",
            "stretch",
            "colormap",
            "remove catalogs",
        ],
    )
    def test_changing_the_display_does_not_resend_the_image(
        self, change, mocker
    ):
        """
        The pixels go to the front end when an image is loaded and at no
        other time. Sending them again for each change of cuts or stretch
        is tens of megabytes a time for a camera frame, which makes a
        slider unusable over a network.
        """
        rng = np.random.default_rng(seed=42)
        arr = rng.integers(1100, 1300, size=(50, 60)).astype(np.uint16)
        self.image.load_image(arr)
        self.image.load_catalog(Table({"x": [1.0], "y": [2.0]}))

        image_mark = self.image._astro_im._image
        sent_before = image_mark.image
        spy_image = mocker.spy(image_mark, "send_state")

        change(self.image)

        # The mark may be asked to send, at the end of a batch, the state
        # that changed in it, but that must never include the image.
        for call in spy_image.call_args_list:
            keys = call.kwargs.get("key", call.args[0] if call.args else None)
            assert keys is not None
            assert "image" not in keys
        assert image_mark.image is sent_before

    @pytest.mark.parametrize(
        "stretch",
        [apviz.LinearStretch(), apviz.SqrtStretch(), apviz.LogStretch()],
        ids=["linear", "sqrt", "log"],
    )
    def test_screen_matches_astropy_scaling(self, stretch):
        """
        What the front end draws, from the pixels, the color limits and
        the colors it is sent, is what astropy's interval and stretch make
        of the image, to within the step between neighboring colors.
        """
        rng = np.random.default_rng(seed=42)
        arr = rng.uniform(1100, 1300, size=(50, 60)).astype(np.float32)
        cuts = apviz.ManualInterval(1120, 1280)

        self.image.load_image(arr)
        self.image.set_cuts(cuts)
        self.image.set_stretch(stretch)
        self.image.set_colormap("gray")

        scaled = cuts(arr)
        expected = stretch(scaled)
        if stretch_fits_in_colors(stretch):
            # A pixel is drawn in the color of the texel it falls in, so
            # it can be off by as much as the stretch changes across that
            # texel.
            texel = 1 / 255
            tolerance = (
                stretch(np.clip(scaled + texel, 0, 1))
                - stretch(np.clip(scaled - texel, 0, 1))
                # The colors themselves are 8 bit.
                + 1 / 255
            )
        else:
            # The stretch is applied to the pixels themselves, so only
            # the rounding to one of 256 colors is left.
            tolerance = 1 / 255
        assert np.all(
            np.abs(screen_levels(self.image) - expected) <= tolerance
        )

    @pytest.mark.parametrize(
        "dtype", [np.uint8, np.uint16, np.int32, np.float32, np.float64]
    )
    def test_image_is_sent_as_float32(self, dtype):
        """
        The front end draws only float32 by value: it reads uint8 as a
        fraction of 255 and refuses the other integer types, so every
        image has to be sent as float32, with its values unchanged.
        """
        arr = np.arange(30 * 40).reshape(30, 40).astype(dtype)

        self.image.load_image(arr)

        sent = self.image._astro_im._image.image
        assert sent.dtype == np.float32
        np.testing.assert_array_equal(sent, arr)

    def test_steep_stretch_is_applied_to_the_pixels(self, mocker):
        """
        A log stretch is too steep at the faint end to be built into 256
        colors: everything in the faintest 1/256 of the cuts would be one
        shade, where the stretch spreads it over a quarter of the way from
        black to white. So the cuts and stretch are applied to the pixels,
        which are sent again, as 8 bit, when either changes.
        """
        rng = np.random.default_rng(seed=42)
        arr = rng.uniform(1100, 1300, size=(50, 60)).astype(np.float32)
        stretch = apviz.LogStretch()
        self.image.load_image(arr)
        image_mark = self.image._astro_im._image
        spy_image = mocker.spy(image_mark, "send_state")

        self.image.set_stretch(stretch)

        assert spy_image.call_count == 1
        assert image_mark.image.dtype == np.uint8
        assert_displayed(self.image, arr, self.image.get_cuts(), stretch)

        cuts = apviz.ManualInterval(1150, 1250)
        self.image.set_cuts(cuts)

        assert spy_image.call_count == 2
        assert_displayed(self.image, arr, cuts, stretch)

    def test_colormap_change_never_resends_the_image(self, mocker):
        """
        A colormap is only ever a list of colors, so changing it must not
        send the image even with a stretch that is applied to the pixels.
        """
        rng = np.random.default_rng(seed=42)
        arr = rng.uniform(1100, 1300, size=(50, 60)).astype(np.float32)
        self.image.load_image(arr)
        self.image.set_stretch(apviz.LogStretch())
        image_mark = self.image._astro_im._image
        sent_before = image_mark.image
        spy_image = mocker.spy(image_mark, "send_state")

        self.image.set_colormap("viridis")

        assert spy_image.call_count == 0
        assert image_mark.image is sent_before
        assert image_mark.scales["image"].colors == bqcolors("viridis")

    def test_back_to_a_gentle_stretch_sends_the_pixels_again(self):
        """
        After a steep stretch the front end holds scaled pixels, so going
        back to a stretch that is built into the colors has to send it
        the pixels as they are once more.
        """
        rng = np.random.default_rng(seed=42)
        arr = rng.uniform(1100, 1300, size=(50, 60)).astype(np.float32)
        self.image.load_image(arr)
        cuts = self.image.get_cuts()
        self.image.set_stretch(apviz.LogStretch())

        self.image.set_stretch(apviz.LinearStretch())

        assert_displayed(self.image, arr, cuts, apviz.LinearStretch())

    def test_float32_image_is_not_copied(self):
        """
        A float32 image is sent as it is. A copy of a camera frame is
        64 MB that a memory-limited server does not have to spare.
        """
        arr = np.arange(30 * 40, dtype=np.float32).reshape(30, 40)

        self.image.load_image(arr)

        assert np.shares_memory(self.image._astro_im._image.image, arr)

    def test_single_valued_image_gets_usable_color_limits(self):
        """
        The front end divides by the difference of the color limits, so
        an image that is one value everywhere must not make them equal.
        It is drawn in the first color, as astropy's intervals scale it
        to zero.
        """
        self.image.load_image(np.full((20, 30), 7.0))
        self.image.set_colormap("gray")

        color_scale = self.image._astro_im._image.scales["image"]
        assert color_scale.min == 7
        assert color_scale.max > color_scale.min
        np.testing.assert_array_equal(screen_levels(self.image), 0)

    def test_cuts_that_are_not_numbers_give_finite_color_limits(self, data):
        """
        The color limits are sent to the front end as JSON, which has no
        NaN, so cuts that come out as NaN must not be sent as they are.
        """
        self.image.load_image(data)

        self.image.set_cuts(apviz.ManualInterval(np.nan, np.nan))

        color_scale = self.image._astro_im._image.scales["image"]
        assert (color_scale.min, color_scale.max) == (0, 1)

    def test_get_viewport_reflects_interactive_pan(self, data):
        # Panning in the browser shifts the bqplot scales directly. Simulate
        # that here and check that get_viewport reports the new center.
        self.image.load_image(data)
        self.image.set_viewport(center=(75, 50), fov=100)

        scales = self.image._astro_im._scales

        # Horizontal pan: shift the x scale. Set min before max so that the
        # center is already correct when the 'max' observer fires.
        dx = 20
        scales["x"].min += dx
        scales["x"].max += dx

        vp = self.image.get_viewport(sky_or_pixel="pixel")
        assert vp["center"][0] == pytest.approx(75 + dx)
        assert vp["center"][1] == pytest.approx(50)

        # Vertical pan: shift the y scale.
        dy = -15
        scales["y"].min += dy
        scales["y"].max += dy

        vp = self.image.get_viewport(sky_or_pixel="pixel")
        assert vp["center"][0] == pytest.approx(75 + dx)
        assert vp["center"][1] == pytest.approx(50 + dy)


@pytest.mark.parametrize(
    "stretch", [apviz.SqrtStretch(), apviz.LogStretch(), apviz.AsinhStretch()]
)
def test_bqcolors_builds_in_a_stretch(stretch):
    """
    With a stretch, color number i of 256 is the color the colormap has at
    the stretch of i / 255: the stretch is applied to how far a pixel is
    between the color limits, which is what the front end looks a color
    up by.
    """
    levels = np.linspace(0, 1, 256)
    greys = [to_rgb(color)[0] for color in bqcolors("gray", stretch=stretch)]

    # The gray colormap rises evenly from black to white, and the colors
    # are 8 bit.
    np.testing.assert_allclose(greys, stretch(levels), atol=1 / 255)
    assert greys[0] == 0
    assert greys[-1] == 1


def test_bqcolors_linear_stretch_changes_nothing():
    """
    A linear stretch leaves the colors exactly as they are with none, so
    that an image nobody has stretched is drawn as it always was.
    """
    assert bqcolors("viridis", stretch=apviz.LinearStretch()) == bqcolors(
        "viridis"
    )


@pytest.mark.parametrize(
    "stretch, fits",
    [
        (apviz.LinearStretch(), True),
        (apviz.SqrtStretch(), True),
        (apviz.AsinhStretch(), True),
        (apviz.LogStretch(), False),
        (apviz.LogStretch(a=10), True),
        (apviz.PowerStretch(0.2), False),
    ],
)
def test_stretch_fits_in_colors(stretch, fits):
    """
    Gentle stretches can be built into 256 colors; one that is steep
    somewhere, like the default log stretch near black, cannot be without
    showing bands.
    """
    assert stretch_fits_in_colors(stretch) is fits


def test_scale_to_levels_is_the_color_each_pixel_is_drawn_in():
    """
    The scaled image is the number, from 0 to 255, of the color astropy's
    interval and stretch put each pixel in, including in an image taller
    than the band of rows that is scaled at a time and one whose values
    run outside the cuts.
    """
    rng = np.random.default_rng(seed=42)
    arr = rng.uniform(900, 1500, size=(600, 37)).astype(np.float32)
    stretch = apviz.LogStretch()

    levels = scale_to_levels(arr, 1000, 1400, stretch)

    assert levels.dtype == np.uint8
    scaled = stretch(apviz.ManualInterval(1000, 1400)(arr))
    np.testing.assert_array_equal(
        levels, np.minimum(np.floor(scaled * 256), 255)
    )
    assert levels.min() == 0
    assert levels.max() == 255


def test_scale_to_levels_keeps_nans():
    """
    The front end leaves NaN pixels undrawn, and 8 bit integers cannot
    hold a NaN, so an image with any is scaled to float32 instead, which
    the front end draws the same way.
    """
    arr = np.linspace(1000, 1400, 600 * 4, dtype=np.float32).reshape(600, 4)
    # In the last band of rows, to check that the bands before it are not
    # what is returned.
    arr[590, 2] = np.nan
    stretch = apviz.LogStretch()

    levels = scale_to_levels(arr, 1000, 1400, stretch)

    assert levels.dtype == np.float32
    assert np.isnan(levels[590, 2])
    np.testing.assert_allclose(
        levels, stretch(apviz.ManualInterval(1000, 1400)(arr))
    )
