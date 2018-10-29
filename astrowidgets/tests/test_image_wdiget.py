import numpy as np
from astropy.table import Table
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
