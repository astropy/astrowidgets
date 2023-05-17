import re

import pytest
import numpy as np
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy.nddata import CCDData

from ..ginga import ImageWidget


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
