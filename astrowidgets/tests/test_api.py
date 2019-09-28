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
    with pytest.raises(AttributeError):
        image.is_marking = True


def test_start_marking():
    image = ImageWidget()

    # Setting these to check that start_marking affects them.
    image.click_center = True
    assert image.click_center
    image.scroll_pan = False
    assert not image.scroll_pan

    marker_style = {'color': 'yellow', 'radius': 10, 'type': 'cross'}
    image.start_marking(marker_name='something',
                        marker=marker_style)
    assert image.is_marking
    assert image.marker == marker_style
    assert not image.click_center
    assert not image.click_drag

    # scroll_pan better activate when marking otherwise there is
    # no way to pan while interactively marking
    assert image.scroll_pan

    # Make sure that when we stop_marking we get our old
    # controls back.
    image.stop_marking()
    assert image.click_center
    assert not image.scroll_pan

    # Make sure that click_drag is restored as expected
    image.click_drag = True
    image.start_marking()
    assert not image.click_drag
    image.stop_marking()
    assert image.click_drag


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
    table = Table(data=np.random.randint(0, 100, [5, 2]),
                  names=['x', 'y'], dtype=('int', 'int'))
    image.add_markers(table, x_colname='x', y_colname='y',
                      skycoord_colname='coord', marker_name='test')
    image.add_markers(table, x_colname='x', y_colname='y',
                      skycoord_colname='coord', marker_name='test2')
    image.reset_markers()
    with pytest.raises(ValueError):
        image.get_markers(marker_name='test')
    with pytest.raises(ValueError):
        image.get_markers(marker_name='test2')


def test_remove_markers():
    image = ImageWidget()
    # Add a tag name...
    image._marktags.add(image._default_mark_tag_name)
    with pytest.raises(ValueError) as e:
        image.remove_markers('arf')
    assert 'arf' in str(e.value)


def test_stretch():
    image = ImageWidget()
    with pytest.raises(ValueError) as e:
        image.stretch = 'not a valid value'
    assert 'must be one of' in str(e.value)

    image.stretch = 'log'
    assert isinstance(image.stretch, (ColorDistBase))


def test_cuts():
    image = ImageWidget()

    # An invalid string should raise an error
    with pytest.raises(ValueError) as e:
        image.cuts = 'not a valid value'
    assert 'must be one of' in str(e.value)

    # Setting cuts to something with incorrect length
    # should raise an error.
    with pytest.raises(ValueError) as e:
        image.cuts = (1, 10, 100)
    assert 'length 2' in str(e.value)

    # These ought to succeed

    image.cuts = 'histogram'
    assert image.cuts == (0.0, 0.0)

    image.cuts = [10, 100]
    assert image.cuts == (10, 100)


def test_colormap():
    image = ImageWidget()
    cmap_desired = 'gray'
    cmap_list = image.colormap_options
    assert len(cmap_list) > 0 and cmap_desired in cmap_list

    image.set_colormap(cmap_desired)


def test_cursor():
    image = ImageWidget()
    assert image.cursor in ALLOWED_CURSOR_LOCATIONS
    with pytest.raises(ValueError):
        image.cursor = 'not a valid option'
    image.cursor = 'bottom'
    assert image.cursor == 'bottom'


def test_click_drag():
    image = ImageWidget()
    # Set this to ensure that click_drag turns it off
    image._click_center = True

    # Make sure that setting click_drag to False does not turn off
    # click_center.

    image.click_drag = False
    assert image.click_center

    image.click_drag = True

    assert not image.click_center

    # If is_marking is true then trying to click_drag
    # should fail.
    image._is_marking = True
    with pytest.raises(ValueError) as e:
        image.click_drag = True
    assert 'Interactive marking' in str(e.value)


def test_click_center():
    image = ImageWidget()
    assert (image.click_center is True) or (image.click_center is False)

    # Set click_drag True and check that click_center affects it appropriately
    image.click_drag = True

    image.click_center = False
    assert image.click_drag

    image.click_center = True
    assert not image.click_drag

    image.start_marking()
    # If marking is in progress then setting click center should fail
    with pytest.raises(ValueError) as e:
        image.click_center = True
    assert 'Cannot set' in str(e.value)

    # setting to False is fine though so no error is expected here
    image.click_center = False


def test_scroll_pan():
    image = ImageWidget()

    # Make sure scroll_pan is actually settable
    for val in [True, False]:
        image.scroll_pan = val
        assert image.scroll_pan is val


def test_save():
    image = ImageWidget()
    filename = 'woot.png'
    image.save(filename)


def test_width_height():
    image = ImageWidget(image_width=250, image_height=100)
    assert image.image_width == 250
    assert image.image_height == 100
