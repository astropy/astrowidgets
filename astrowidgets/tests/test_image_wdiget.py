import numpy as np
from astropy.table import Table
from astropy.wcs import WCS
from astropy.nddata import CCDData
from astropy.coordinates import SkyCoord

from ..core import ImageWidget


def test_setting_image_width_height():
    image = ImageWidget()
    width = 200
    height = 300
    image.image_width = width
    image.image_height = height
    assert image._viewer.get_window_size() == (width, height)


def test_add_marker_does_not_modify_input_table():
    # Regression test for #45
    # Adding markers should not modify the input data table
    image = ImageWidget(image_width=300, image_height=300)
    data = np.random.random([300, 300])
    image.load_array(data)
    x = [20, 30, 40]
    y = [40, 80, 100]
    # Create two separate tables for comparison after add_markers.
    orig_table = Table(data=[x, y], names=['x', 'y'])
    in_table = Table(data=[x, y], names=['x', 'y'])
    image.add_markers(in_table, pixel_coords_offset=5)
    assert (in_table == orig_table).all()


def test_adding_markers_as_world_recovers_with_get_markers():
    """
    Make sure that our internal conversion from world to pixel
    coordinates doesn't mess anything up.
    """
    fake_image = np.random.randn(2000, 2000)
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = (fake_image.shape[0] / 2, fake_image.shape[1] / 2)
    wcs.wcs.ctype = ('RA---TAN', 'DEC--TAN')
    wcs.wcs.crval = (314.275419158, 31.6662781301)
    wcs.wcs.pc = [[0.000153051015113, -3.20700931602e-05],
                  [3.20704370872e-05, 0.000153072382405]]
    fake_ccd = CCDData(data=fake_image, wcs=wcs, unit='adu')
    iw = ImageWidget()
    iw.load_nddata(fake_ccd)
    # Get me 1000 positions please
    marker_locs = np.random.randint(10, high=1990, size=(1000, 2))
    marks_pix = Table(data=marker_locs, names=['x', 'y'])
    marks_world = wcs.all_pix2world(marker_locs, 0)
    marks_coords = SkyCoord(marks_world, unit='degree')
    mark_coord_table = Table(data=[marks_coords], names=['coord'])
    iw.add_markers(mark_coord_table, use_skycoord=True)
    result = iw.get_markers(pixel_coords_offset=0)
    # Check the x, y positions as long as we are testing things...
    np.testing.assert_allclose(result['x'], marks_pix['x'])
    np.testing.assert_allclose(result['y'], marks_pix['y'])
    np.testing.assert_allclose(result['coord'].ra.deg,
                               mark_coord_table['coord'].ra.deg)
    np.testing.assert_allclose(result['coord'].dec.deg,
                               mark_coord_table['coord'].dec.deg)
