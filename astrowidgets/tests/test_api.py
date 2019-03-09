

import numpy as np

import pytest

from astropy.io import fits
from astropy.nddata import NDData
from astropy.table import Table

from ginga.ColorDist import ColorDistBase

from ..core import ImageWidget, ALLOWED_CURSOR_LOCATIONS


def test_load_fits():
    image = ImageWidget()
    data = np.random.random([100, 100])
    hdu = fits.PrimaryHDU(data=data)
    image.load_fits(hdu)


def test_load_nddata():
    image = ImageWidget()
    data = np.random.random([100, 100])
    nddata = NDData(data)
    image.load_nddata(nddata)


def test_load_array():
    image = ImageWidget()
    data = np.random.random([100, 100])
    image.load_array(data)


def test_center_on():
    image = ImageWidget()
    x = 10
    y = 10
    image.center_on((x, y))


def test_offset_to():
    image = ImageWidget()
    dx = 10
    dy = 10
    image.offset_to(dx, dy)


def test_zoom_level():
    image = ImageWidget()
    image.zoom_level = 5
    assert image.zoom_level == 5


def test_zoom():
    image = ImageWidget()
    image.zoom_level = 3
    val = 2
    image.zoom(val)
    assert image.zoom_level == 6


@pytest.mark.xfail(reason='Not implemented yet')
def test_select_points():
    image = ImageWidget()
    image.select_points()


def test_get_selection():
    image = ImageWidget()
    marks = image.get_markers()
    assert isinstance(marks, Table) or marks is None


def test_stop_marking():
    image = ImageWidget()
    # This is not much of a test...
    image.stop_marking(clear_markers=True)
    assert image.get_markers() is None
    assert image.is_marking is False


def test_is_marking():
    image = ImageWidget()
    assert image.is_marking in [True, False]


def test_add_markers():
    image = ImageWidget()
    table = Table(data=np.random.randint(0, 100, [5, 2]),
                  names=['x', 'y'], dtype=('int', 'int'))
    image.add_markers(table, x_colname='x', y_colname='y',
                      skycoord_colname='coord')


def test_set_markers():
    image = ImageWidget()
    image.marker = {'color': 'yellow', 'radius': 10, 'type': 'cross'}
    assert 'cross' in str(image.marker)
    assert 'yellow' in str(image.marker)
    assert '10' in str(image.marker)


def test_reset_markers():
    image = ImageWidget()
    # First test: this shouldn't raise any errors
    # (it also doesn't *do* anything...)
    image.reset_markers()
    assert image.get_markers() is None

    # TODO add test of actually removing markers...


def test_remove_markers():
    image = ImageWidget()
    # Add a tag name...
    image._marktags.add(image._default_mark_tag_name)
    with pytest.raises(ValueError) as e:
        image.remove_markers('arf')
    assert 'arf' in str(e)


def test_stretch():
    image = ImageWidget()
    with pytest.raises(ValueError) as e:
        image.stretch = 'not a valid value'
    assert 'must be one of' in str(e)

    image.stretch = 'log'
    assert isinstance(image.stretch, (ColorDistBase))


def test_cuts():
    image = ImageWidget()

    # An invalid string should raise an error
    with pytest.raises(ValueError) as e:
        image.cuts = 'not a valid value'
    assert 'must be one of' in str(e)

    # Setting cuts to something with incorrect length
    # should raise an error.
    with pytest.raises(ValueError) as e:
        image.cuts = (1, 10, 100)
    assert 'length 2' in str(e)

    # These ought to succeed

    image.cuts = 'histogram'
    assert image.cuts == (0.0, 0.0)

    image.cuts = [10, 100]
    assert image.cuts == (10, 100)


def test_cursor():
    image = ImageWidget()
    assert image.cursor in ALLOWED_CURSOR_LOCATIONS
    with pytest.raises(ValueError):
        image.cursor = 'not a valid option'
    image.cursor = 'bottom'
    assert image.cursor == 'bottom'


def test_click_drag():
    image = ImageWidget()
    with pytest.raises(NotImplementedError):
        image.click_drag()


def test_click_center():
    image = ImageWidget()
    assert (image.click_center is True) or (image.click_center is False)


def test_scroll_pan():
    image = ImageWidget()
    with pytest.raises(NotImplementedError):
        image.scroll_pan()


def test_save():
    image = ImageWidget()
    filename = 'woot.png'
    image.save(filename)


def test_width_height():
    image = ImageWidget(image_width=250, image_height=100)
    assert image.image_width == 250
    assert image.image_height == 100
