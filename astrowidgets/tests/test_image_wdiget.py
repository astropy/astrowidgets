from ..core import ImageWidget


def test_setting_image_width_height():
    image = ImageWidget()
    width = 200
    height = 300
    image.image_width = width
    image.image_height = height
    assert image._viewer.get_window_size() == (width, height)
