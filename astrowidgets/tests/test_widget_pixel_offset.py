# TODO: How to enable switching out backend and still run the same tests?

import pytest

ginga = pytest.importorskip("ginga")

import re  # noqa: E402

import numpy as np  # noqa: E402

from ginga.misc.log import get_logger  # noqa: E402

from astrowidgets.ginga import ImageWidget  # noqa: E402


class TestGingaWidgetWithPixelOffset:
    def setup_class(self):
        # The pixel offset below is nonsensical. It is chosen simply
        # to make it easy to check for.
        self.offset = 3

        logger = get_logger('my_viewer', log_stderr=False, level=30)
        self.image = ImageWidget(logger, image_width=300, image_height=300,
                                 pixel_coords_offset=self.offset)

        rng = np.random.default_rng(1234)
        data = rng.random((300, 300))
        self.image.load_array(data)

    def test_offset_value(self):
        # Ensure it cannot change after init
        with pytest.raises(AttributeError):
            self.image.pixel_offset = 0

        assert self.image.pixel_offset == self.offset

    def test_move_callback_includes_offset(self):
        # Send a fake move to the callback. What gets put in the cursor
        # value should be the event we sent in plus the offset.
        self.image.click_center = True
        data_x = 100
        data_y = 200
        self.image._mouse_move_cb(self.image.viewer, None, data_x, data_y)
        output_contents = self.image._jup_coord.value
        x_out = re.search(r'X: ([\d\.\d]+)', output_contents)
        x_out = x_out.groups(1)[0]
        y_out = re.search(r'Y: ([\d\.\d]+)', output_contents)
        y_out = y_out.groups(1)[0]
        assert float(x_out) == data_x + self.offset
        assert float(y_out) == data_y + self.offset
