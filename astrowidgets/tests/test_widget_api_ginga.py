import warnings

import numpy as np
import pytest

import ipywidgets as ipyw
from astropy.coordinates import SkyCoord
from astropy.nddata import NDData
from astropy.utils.exceptions import AstropyUserWarning
from astropy.visualization import (AsinhStretch, ContrastBiasStretch,
                                   LinearStretch, LogStretch, PowerDistStretch,
                                   PowerStretch, SinhStretch, SqrtStretch,
                                   SquaredStretch)
from astropy.wcs import WCS
from traitlets import TraitError

from astro_image_display_api.api_test import ImageAPITest
from astro_image_display_api import ImageViewerInterface


def _make_wcs():
    w = WCS(naxis=2)
    w.wcs.crpix = [-234.75, 8.3393]
    w.wcs.cdelt = np.array([-0.066667, 0.066667])
    w.wcs.crval = [0, -90]
    w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    w.wcs.set_pv([(2, 1, 45.0)])
    return w


_ = pytest.importorskip("ginga",
                        reason="Package required for test is not "
                               "available.")
from ginga import ColorDist  # noqa: E402
from astrowidgets.ginga import ImageWidget  # noqa: E402


def _loaded_widget():
    image = ImageWidget()
    image.load_image(np.random.default_rng(1234).random((100, 150)))
    return image


def _astropy_levels(stretch, dist):
    # Reproduce the integer lookup table ginga would build from an astropy
    # stretch evaluated on the same 0..1 ramp, so it can be compared against a
    # configured ginga ColorDist's ``hash``.
    base = np.arange(0.0, float(dist.hashsize), 1.0) / dist.hashsize
    out = np.clip(np.asarray(stretch(base, clip=True), dtype=float), 0.0, 1.0)
    return (out * (dist.colorlen - 1)).astype(np.int64)


def _max_level_diff(dist, stretch):
    return int(np.abs(dist.hash.astype(np.int64)
                      - _astropy_levels(stretch, dist)).max())


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


def test_pixel_offset_removed():
    # The pixel_coords_offset kwarg and pixel_offset property have been
    # removed; the API-approved way to handle this is via the viewport.
    image = ImageWidget()
    assert not hasattr(image, 'pixel_offset')
    # pixel_coords_offset is no longer a recognized parameter, so it falls
    # through to traitlets as an unrecognized argument.
    with pytest.warns(DeprecationWarning, match='pixel_coords_offset'):
        ImageWidget(pixel_coords_offset=1)


def test_image_size_getters_are_read_only():
    # image_width/image_height are read-only; resizing goes through the
    # viewport API instead.
    image = ImageWidget(image_width=400, image_height=300)
    assert isinstance(image.image_width, int)
    assert isinstance(image.image_height, int)
    assert image.image_width == 400
    assert image.image_height == 300
    with pytest.raises(AttributeError):
        image.image_width = 600
    with pytest.raises(AttributeError):
        image.image_height = 600


@pytest.mark.parametrize('stretch', [AsinhStretch(0.05), SinhStretch(0.2),
                                     SqrtStretch(), SquaredStretch(),
                                     LinearStretch()])
def test_native_dist_matches_astropy_exactly(stretch):
    # For the families ginga implements directly, the configured ginga
    # ColorDist is an exact reparametrization of the astropy stretch, so its
    # lookup table is identical to evaluating the astropy stretch on the ramp.
    image = _loaded_widget()
    image.set_stretch(stretch)
    dist = image._viewer.get_rgbmap().get_dist()
    np.testing.assert_array_equal(dist.hash, _astropy_levels(stretch, dist))


@pytest.mark.parametrize('stretch,max_diff',
                         [(LogStretch(500), 1), (PowerDistStretch(200), 2)])
def test_native_dist_matches_astropy_closely(stretch, max_diff):
    # Ginga normalizes log/power slightly differently than astropy (log(a) vs
    # log(a+1)); the curve shape still matches to within a quantization level.
    image = _loaded_widget()
    image.set_stretch(stretch)
    dist = image._viewer.get_rgbmap().get_dist()
    assert _max_level_diff(dist, stretch) <= max_diff


def test_asinh_parameter_flows_to_ginga():
    # The astropy shape parameter must reach ginga, not be discarded:
    # AsinhStretch(a) -> AsinhDist(factor=1/a, nonlinearity=arcsinh(1/a)).
    image = _loaded_widget()
    image.set_stretch(AsinhStretch(0.05))
    dist = image._viewer.get_rgbmap().get_dist()
    assert isinstance(dist, ColorDist.AsinhDist)
    assert dist.factor == pytest.approx(20.0)
    assert dist.nonlinearity == pytest.approx(np.arcsinh(20.0))


def test_log_parameter_flows_to_ginga():
    image = _loaded_widget()
    image.set_stretch(LogStretch(500))
    dist = image._viewer.get_rgbmap().get_dist()
    assert isinstance(dist, ColorDist.LogDist)
    assert dist.exp == 500


def test_power_stretch_uses_adapter_not_ginga_power():
    # ginga's 'power' algorithm is (a**x - 1)/(a - 1) -- that is astropy's
    # PowerDistStretch, NOT PowerStretch (x**a). PowerStretch has no ginga
    # family and must be applied faithfully through the adapter.
    image = _loaded_widget()
    stretch = PowerStretch(3.0)
    image.set_stretch(stretch)
    dist = image._viewer.get_rgbmap().get_dist()
    assert not isinstance(dist, ColorDist.PowerDist)
    np.testing.assert_array_equal(dist.hash, _astropy_levels(stretch, dist))


def test_powerdist_stretch_maps_to_ginga_power():
    image = _loaded_widget()
    image.set_stretch(PowerDistStretch(200))
    dist = image._viewer.get_rgbmap().get_dist()
    assert isinstance(dist, ColorDist.PowerDist)
    assert dist.exp == 200


def test_ginga_less_stretch_applied_faithfully_without_warning():
    # ContrastBiasStretch has no ginga family; the adapter applies it exactly
    # and, unlike the old behavior, does not warn or fall back to linear.
    image = _loaded_widget()
    stretch = ContrastBiasStretch(1.0, 0.5)
    with warnings.catch_warnings():
        warnings.simplefilter('error', AstropyUserWarning)
        image.set_stretch(stretch)
    dist = image._viewer.get_rgbmap().get_dist()
    np.testing.assert_array_equal(dist.hash, _astropy_levels(stretch, dist))


@pytest.mark.parametrize('stretch,algorithm',
                         [(LinearStretch(), 'linear'), (LogStretch(), 'log'),
                          (AsinhStretch(), 'asinh'), (SinhStretch(), 'sinh'),
                          (SqrtStretch(), 'sqrt'), (SquaredStretch(), 'squared')])
def test_native_stretch_reports_family_name(stretch, algorithm):
    # Native families keep ginga's algorithm-name introspection working and
    # apply without warning.
    image = _loaded_widget()
    with warnings.catch_warnings():
        warnings.simplefilter('error', AstropyUserWarning)
        image.set_stretch(stretch)
    assert image._viewer.get_rgbmap().get_hash_algorithm() == algorithm


@pytest.mark.parametrize('fov', [50, 120, 200])
def test_non_square_viewport_fov_round_trips(fov):
    # The widget need not be square: fov is defined as the horizontal field of
    # view, so it round-trips through set/get_viewport for a non-square viewer.
    image = ImageWidget(image_width=600, image_height=300)
    image.load_image(np.random.default_rng(1234).random((100, 150)))
    image.set_viewport(center=(40, 60), fov=fov)
    vport = image.get_viewport(sky_or_pixel='pixel')
    assert vport['fov'] == pytest.approx(fov)
    assert vport['center'][0] == pytest.approx(40)
    assert vport['center'][1] == pytest.approx(60)


def test_center_on_pixel_updates_viewport():
    # _center_on with a pixel tuple should re-center the stored viewport.
    image = ImageWidget()
    image.load_image(np.zeros((100, 150)))
    image._center_on((30, 40))
    center = image.get_viewport(sky_or_pixel='pixel')['center']
    assert center[0] == pytest.approx(30)
    assert center[1] == pytest.approx(40)


def test_center_on_skycoord_updates_viewport():
    # _center_on with a SkyCoord should re-center the stored viewport.
    wcs = _make_wcs()
    image = ImageWidget()
    image.load_image(NDData(data=np.zeros((100, 150)), wcs=wcs))
    target = SkyCoord(*wcs.wcs.crval, unit='deg')
    image._center_on(target)
    center = image.get_viewport(sky_or_pixel='sky')['center']
    assert isinstance(center, SkyCoord)
    assert center.separation(target).arcsec == pytest.approx(0, abs=1e-3)


def _two_image_widget():
    # A widget with two labeled images; 'b' is the one currently displayed.
    rng = np.random.default_rng(1234)
    image = ImageWidget()
    image.load_image(rng.random((100, 100)), image_label='a')
    image.load_image(rng.random((100, 100)), image_label='b')
    return image


def test_get_viewport_does_not_corrupt_non_displayed_state():
    # Reading the viewport of a non-displayed image must not "sync" it from
    # the live ginga pan/scale, which belongs to a different image.
    rng = np.random.default_rng(1234)
    image = ImageWidget()
    image.load_image(rng.random((256, 256)), image_label='a')
    image.set_viewport(center=(10, 20), fov=50, image_label='a')
    # Load a differently-sized image so the live viewer state clearly does
    # not match the stored viewport for 'a'.
    image.load_image(rng.random((300, 200)), image_label='b')

    # Read twice: the first call must return the stored values, and it must
    # not overwrite them, so the second call agrees.
    for _ in range(2):
        vport = image.get_viewport(sky_or_pixel='pixel', image_label='a')
        assert vport['center'][0] == pytest.approx(10)
        assert vport['center'][1] == pytest.approx(20)
        assert vport['fov'] == pytest.approx(50)


def test_set_cuts_non_displayed_label_leaves_display_alone():
    image = _two_image_widget()

    displayed_before = image._viewer.get_cut_levels()
    # 'a' is not displayed, so the live cut levels must not change...
    image.set_cuts((0, 0.5), image_label='a')
    assert image._viewer.get_cut_levels() == displayed_before
    # ...but the stored cuts for 'a' are still updated.
    stored = image.get_cuts(image_label='a')
    assert (stored.vmin, stored.vmax) == (0, 0.5)

    # The same call for the displayed image does take effect.
    image.set_cuts((0, 0.5), image_label='b')
    assert image._viewer.get_cut_levels() == pytest.approx((0, 0.5))


def test_set_stretch_non_displayed_label_leaves_display_alone():
    image = _two_image_widget()

    dist_before = image._viewer.get_rgbmap().get_dist()
    # 'a' is not displayed, so the live color distribution must not change.
    image.set_stretch(LogStretch(500), image_label='a')
    assert image._viewer.get_rgbmap().get_dist() is dist_before

    # The same call for the displayed image does take effect.
    image.set_stretch(LogStretch(500), image_label='b')
    dist = image._viewer.get_rgbmap().get_dist()
    assert isinstance(dist, ColorDist.LogDist)
    assert dist.exp == 500


def test_set_colormap_non_displayed_label_leaves_display_alone():
    image = _two_image_widget()

    cmap_before = image._viewer.get_settings().get_setting('color_map').value
    assert cmap_before != 'viridis'
    # 'a' is not displayed, so the live colormap must not change...
    image.set_colormap('viridis', image_label='a')
    assert (image._viewer.get_settings().get_setting('color_map').value
            == cmap_before)
    # ...but the stored colormap for 'a' is still updated.
    assert image.get_colormap(image_label='a') == 'viridis'

    # The same call for the displayed image does take effect.
    image.set_colormap('viridis', image_label='b')
    assert (image._viewer.get_settings().get_setting('color_map').value
            == 'viridis')


def test_set_viewport_non_displayed_label_leaves_display_alone():
    image = _two_image_widget()

    pan_before = image._viewer.get_pan()
    scale_before = image._viewer.get_scale()
    # 'a' is not displayed, so the live pan/scale must not change...
    image.set_viewport(center=(5, 5), fov=10, image_label='a')
    assert image._viewer.get_pan() == pytest.approx(pan_before)
    assert image._viewer.get_scale() == pytest.approx(scale_before)
    # ...but the stored viewport for 'a' is still updated.
    vport = image.get_viewport(sky_or_pixel='pixel', image_label='a')
    assert vport['center'][0] == pytest.approx(5)
    assert vport['center'][1] == pytest.approx(5)
    assert vport['fov'] == pytest.approx(10)


def test_set_colormap_invalid_name_raises():
    image = _loaded_widget()
    with pytest.raises(ValueError, match='not a valid ginga'):
        image.set_colormap('not-a-real-colormap')


class TestGingaWidget(ImageAPITest):
    image_widget_class = ImageWidget
    cursor_error_classes = (ValueError, TraitError)
