import numpy as np
import pytest

import ipywidgets as ipyw
from astropy.coordinates import SkyCoord
from astropy.nddata import NDData
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


class TestGingaWidget(ImageAPITest):
    image_widget_class = ImageWidget
    cursor_error_classes = (ValueError, TraitError)
