import re

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
    image = ImageWidget(image_width=300, image_height=300,
                        pixel_coords_offset=5)
    data = np.random.random([300, 300])
    image.load_array(data)
    x = [20, 30, 40]
    y = [40, 80, 100]
    # Create two separate tables for comparison after add_markers.
    orig_table = Table(data=[x, y], names=['x', 'y'])
    in_table = Table(data=[x, y], names=['x', 'y'])
    image.add_markers(in_table)
    assert (in_table == orig_table).all()


def test_adding_markers_as_world_recovers_with_get_markers():
    """
    Make sure that our internal conversion from world to pixel
    coordinates doesn't mess anything up.
    """
    npix_side = 100
    fake_image = np.random.randn(npix_side, npix_side)
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = (fake_image.shape[0] / 2, fake_image.shape[1] / 2)
    wcs.wcs.ctype = ('RA---TAN', 'DEC--TAN')
    wcs.wcs.crval = (314.275419158, 31.6662781301)
    wcs.wcs.pc = [[0.000153051015113, -3.20700931602e-05],
                  [3.20704370872e-05, 0.000153072382405]]
    fake_ccd = CCDData(data=fake_image, wcs=wcs, unit='adu')
    iw = ImageWidget(pixel_coords_offset=0)
    iw.load_nddata(fake_ccd)
    # Get me 100 positions please, not right at the edge
    marker_locs = np.random.randint(10,
                                    high=npix_side - 10,
                                    size=(100, 2))
    marks_pix = Table(data=marker_locs, names=['x', 'y'])
    marks_world = wcs.all_pix2world(marker_locs, 0)
    marks_coords = SkyCoord(marks_world, unit='degree')
    mark_coord_table = Table(data=[marks_coords], names=['coord'])
    iw.add_markers(mark_coord_table, use_skycoord=True)
    result = iw.get_markers()
    # Check the x, y positions as long as we are testing things...
    np.testing.assert_allclose(result['x'], marks_pix['x'])
    np.testing.assert_allclose(result['y'], marks_pix['y'])
    np.testing.assert_allclose(result['coord'].ra.deg,
                               mark_coord_table['coord'].ra.deg)
    np.testing.assert_allclose(result['coord'].dec.deg,
                               mark_coord_table['coord'].dec.deg)


def test_can_set_pixel_offset_at_object_level():
    # The pixel offset below is nonsensical. It is chosen simply
    # to make it easy to check for.
    offset = 3
    image = ImageWidget(image_width=300, image_height=300,
                        pixel_coords_offset=offset)
    assert image._pixel_offset == offset


def test_move_callback_includes_offset():
    # The pixel offset below is nonsensical. It is chosen simply
    # to make it easy to check for.
    offset = 3
    image = ImageWidget(image_width=300, image_height=300,
                        pixel_coords_offset=offset)
    data = np.random.random([300, 300])
    image.load_array(data)
    # Send a fake move to the callback. What gets put in the cursor
    # value should be the event we sent in plus the offset.
    image.click_center = True
    data_x = 100
    data_y = 200
    image._mouse_move_cb(image._viewer, None, data_x, data_y)
    output_contents = image._jup_coord.value
    print(output_contents)
    x_out = re.search(r'X: ([\d\.\d]+)', output_contents)
    x_out = x_out.groups(1)[0]
    y_out = re.search(r'Y: ([\d\.\d]+)', output_contents)
    y_out = y_out.groups(1)[0]
    assert float(x_out) == data_x + offset
    assert float(y_out) == data_y + offset
    # image.print_out.get_state()['outputs']
