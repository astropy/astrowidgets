import re

import pytest
import numpy as np
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy.nddata import CCDData
from astropy.coordinates import SkyCoord

from ..core import ImageWidget, RESERVED_MARKER_SET_NAMES


def _make_fake_ccd(with_wcs=True):
    """
    Generate a CCDData object for use with ImageWidget tests.

    Parameters
    ----------

    with_wcs : bool, optional
        If ``True`` the image will have a WCS attached to it,
        which is useful for some of the marker tests.

    Returns
    -------

    `astropy.nddata.CCDData`
        CCD image
    """
    npix_side = 100
    fake_image = np.random.randn(npix_side, npix_side)
    if with_wcs:
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = (fake_image.shape[0] / 2, fake_image.shape[1] / 2)
        wcs.wcs.ctype = ('RA---TAN', 'DEC--TAN')
        wcs.wcs.crval = (314.275419158, 31.6662781301)
        wcs.wcs.pc = [[0.000153051015113, -3.20700931602e-05],
                      [3.20704370872e-05, 0.000153072382405]]
    else:
        wcs = None

    return CCDData(data=fake_image, wcs=wcs, unit='adu')


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
    fake_ccd = _make_fake_ccd(with_wcs=True)
    npix_side = fake_ccd.shape[0]
    wcs = fake_ccd.wcs
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
    x_out = re.search(r'X: ([\d\.\d]+)', output_contents)
    x_out = x_out.groups(1)[0]
    y_out = re.search(r'Y: ([\d\.\d]+)', output_contents)
    y_out = y_out.groups(1)[0]
    assert float(x_out) == data_x + offset
    assert float(y_out) == data_y + offset


def test_can_add_markers_with_names():
    """
    Test a few things related to naming marker sets
    """
    npix_side = 200
    image = ImageWidget(image_width=npix_side,
                        image_height=npix_side)
    x = np.array([20, 30, 40])
    y = np.array([40, 80, 100])

    # This should succeed without error
    image.add_markers(Table(data=[x, y], names=['x', 'y']),
                      marker_name='nonsense')

    # The name 'nonsense', and nothing else, should be in the
    # set of markers.
    assert set(['nonsense']) == image._marktags

    # Add more markers with the same name
    # This should succeed without error
    image.add_markers(Table(data=[x, y], names=['x', 'y']),
                      marker_name='nonsense')

    # check that we get the right number of markers
    marks = image.get_markers(marker_name='nonsense')
    assert len(marks) == 6

    # Make sure setting didn't change the default name
    assert image._default_mark_tag_name == 'default-marker-name'

    # Try adding markers without a name
    image.add_markers(Table(data=[x, y], names=['x', 'y']))
    assert image._marktags == set(['nonsense', image._default_mark_tag_name])

    # Delete just the nonsense markers
    image.remove_markers('nonsense')

    assert 'nonsense' not in image._marktags
    assert image._default_mark_tag_name in image._marktags

    # Add the nonsense markers back...
    image.add_markers(Table(data=[x, y], names=['x', 'y']),
                      marker_name='nonsense')
    # ...and now delete all of the markers
    image.reset_markers()
    # We should have no markers on the image
    assert image._marktags == set()

    # Simulate a mouse click and make sure the expected marker
    # name has been added.
    data_x = 50
    data_y = 50
    image._is_marking = True
    image._mouse_click_cb(image._viewer, None, data_x, data_y)
    assert image._interactive_marker_set_name in image._marktags


def test_mark_with_reserved_name_raises_error():
    npix_side = 200
    image = ImageWidget(image_width=npix_side,
                        image_height=npix_side)
    x = np.array([20, 30, 40])
    y = np.array([40, 80, 100])
    for name in RESERVED_MARKER_SET_NAMES:
        with pytest.raises(ValueError):
            image.add_markers(Table(data=[x, y], names=['x', 'y']),
                              marker_name=name)


def test_get_marker_with_names():
    # Check a few ways of getting markers out
    npix_side = 200
    image = ImageWidget(image_width=npix_side,
                        image_height=npix_side)

    x = np.array([20, 30, 40])
    y = np.array([40, 80, 100])
    input_markers = Table(data=[x, y], names=['x', 'y'])
    # Add some markers with our own name
    image.add_markers(input_markers, marker_name='nonsense')

    # Add same markers without a name so that name defaults to
    # image._default_mark_tag_name
    image.add_markers(input_markers)

    # Add pseudo-interactive points
    image._is_marking = True
    for data_x, data_y in input_markers:
        image._mouse_click_cb(image._viewer, None, data_x, data_y)

    # Should have three sets of markers: nonsense, default non-interactive,
    # interactive
    assert len(image._marktags) == 3

    for marker in image._marktags:
        out_table = image.get_markers(marker_name=marker)
        # No guarantee markers will come back in the same order, so sort them.
        out_table.sort('x')
        assert (out_table['x'] == input_markers['x']).all()
        assert (out_table['y'] == input_markers['y']).all()

    # Get all of markers at once
    all_marks = image.get_markers(marker_name='all')

    # That should have given us three copies of the input table
    expected = vstack([input_markers] * 3, join_type='exact')

    # Sort before comparing
    expected.sort(['x', 'y'])
    all_marks.sort(['x', 'y'])

    assert (expected['x'] == all_marks['x']).all()
    assert (expected['y'] == all_marks['y']).all()


def test_unknown_marker_name_error():
    """
    Regression test for https://github.com/astropy/astrowidgets/issues/97

    This particular test checks that getting a marker name that
    does not exist raises an error.
    """
    iw = ImageWidget()
    bad_name = 'not a real marker name'
    with pytest.raises(ValueError) as e:
        iw.get_markers(marker_name=bad_name)

    assert f"No markers named '{bad_name}'" in str(e.value)


def test_marker_name_has_no_marks_warning():
    """
    Regression test for https://github.com/astropy/astrowidgets/issues/97

    This particular test checks that getting an empty table gives a
    useful warning message.
    """
    iw = ImageWidget()
    bad_name = 'empty marker set'
    iw.start_marking(marker_name=bad_name)

    with pytest.warns(UserWarning) as record:
        iw.get_markers(marker_name=bad_name)

    assert f"Marker set named '{bad_name}' is empty" in str(record[0].message)


def test_empty_marker_name_works_with_all():
    """
    Regression test for https://github.com/astropy/astrowidgets/issues/97

    This particular test checks that an empty table doesn't break
    marker_name='all'. The bug only comes up if there is a coordinate
    column, so use a fake image a WCS.
    """
    iw = ImageWidget()
    fake_ccd = _make_fake_ccd(with_wcs=True)
    iw.load_nddata(fake_ccd)

    x = np.array([20, 30, 40])
    y = np.array([40, 80, 100])
    input_markers = Table(data=[x, y], names=['x', 'y'])
    # Add some markers with our own name
    iw.add_markers(input_markers, marker_name='nonsense')

    # Start marking to create a new marker set that is empty
    iw.start_marking(marker_name='empty')

    marks = iw.get_markers(marker_name='all')
    assert len(marks) == len(x)
    assert 'empty' not in marks['marker name']


def test_add_single_marker():
    """
    Test a few things related to naming marker sets
    """
    fake_ccd = _make_fake_ccd(with_wcs=True)
    npix_side = fake_ccd.shape[0]
    wcs = fake_ccd.wcs
    iw = ImageWidget(pixel_coords_offset=0)
    iw.load_nddata(fake_ccd)
    # Get me 100 positions please, not right at the edge
    marker_locs = np.random.randint(10,
                                    high=npix_side - 10,
                                    size=(100, 2))
    marks_world = wcs.all_pix2world(marker_locs, 0)
    marks_coords = SkyCoord(marks_world, unit='degree')
    mark_coord_table = Table(data=[marks_coords], names=['coord'])
    iw.add_markers(mark_coord_table[0], use_skycoord=True)
