import numpy as np
from astropy.io import fits
from astropy.nddata import NDData
from astropy.table import Table

from ..core import ImageWidget


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
    image.center_on(x, y)


def test_offset_to():
    image = ImageWidget()
    dx = 10
    dy = 10
    image.offset_to(dx, dy)


def test_zoom_level():
    image = ImageWidget()
    image.zoom_level()


def test_zoom():
    image = ImageWidget()
    val = 2
    image.zoom(val)


def test_select_points():
    image = ImageWidget()
    image.select_points()


def test_get_selection():
    image = ImageWidget()
    image.get_selection()


def test_stop_selecting():
    image = ImageWidget()
    image.stop_selecting(clear_marks=True)


def test_is_selecting():
    image = ImageWidget()
    image.is_selecting()


def test_add_marks():
    image = ImageWidget()
    table = Table(data=np.random.randint(0, 100, [2, 5]),
                  names=['x', 'y'])
    image.add_marks(table, x_colname='x', y_colname='y',
                    skycoord_colname='coord')


def test_reset_marks():
    image = ImageWidget()
    image.reset_marks()


def test_stretch():
    image = ImageWidget()
    image.stretch()


def test_cuts():
    image = ImageWidget()
    image.cuts()


def test_cursor():
    image = ImageWidget()
    image.cursor()


def test_click_drag():
    image = ImageWidget()
    image.click_drag()


def test_click_center():
    image = ImageWidget()
    image.click_center()


def test_scroll_pan():
    image = ImageWidget()
    image.scroll_pan()


def test_save():
    image = ImageWidget()
    filename = 'woot.png'
    image.save(filename)
