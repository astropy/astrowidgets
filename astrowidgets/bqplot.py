"""``astrowidgets`` with ``bqplot`` as backend."""

import numpy as np

# Jupyter widgets
import ipywidgets as ipyw

# bqplot
from bqplot import Figure, HeatMap, LinearScale, ColorScale, Axis, ColorAxis

__all__ = ['ImageWidget']


class ImageWidget(ipyw.VBox):
    """Image widget for Jupyter notebook using ``bqplot``.

    .. todo:: Any property passed to constructor has to be valid keyword.

    Parameters
    ----------
    logger : obj or ``None``
        Python logger.

    image_width, image_height : int
        Dimension of Jupyter notebook's image widget.

    pixel_coords_offset : int, optional
        An offset, typically either 0 or 1, to add/subtract to all
        pixel values when going to/from the displayed image.
        *In almost all situations the default value, ``0``, is the
        correct value to use.*

    """

    def __init__(self, logger=None, image_width=500, image_height=500,
                 pixel_coords_offset=0):
        super().__init__()

        self._logger = logger
        self._pixel_offset = pixel_coords_offset
        self._jup_img = Figure(layout=ipyw.Layout(width=f'{image_width}px',
                                                  height=f'{image_height}px'),
                               padding_y=0)

        # coordinates display
        self._jup_coord = ipyw.HTML('Coordinates show up here')

        self.children = [self._jup_img, self._jup_coord]

    @property
    def logger(self):
        """Logger for this widget."""
        return self._logger

    @property
    def image_width(self):
        return int(self._jup_img.layout.width.replace('px', ''))

    @image_width.setter
    def image_width(self, value):
        # widgets expect width/height as strings, but most users will not, so
        # do the conversion.
        self._jup_img.layout.width = f'{value}px'

    @property
    def image_height(self):
        return int(self._jup_img.layout.height.replace('px', ''))

    @image_height.setter
    def image_height(self, value):
        # widgets expect width/height as strings, but most users will not, so
        # do the conversion.
        self._jup_img.layout.height = f'{value}px'

    @property
    def pixel_offset(self):
        """An offset, typically either 0 or 1, to add/subtract to all
        pixel values when going to/from the displayed image.
        *In almost all situations the default value, ``0``, is the
        correct value to use.*

        This value cannot be modified after initialization.
        """
        return self._pixel_offset

    def _mouse_move_cb(self, *args):
        """
        Callback to display position in RA/DEC deg.
        """
        # FIXME: Where is the X and Y info?
        self._jup_coord.value = f'{args}'

        # TODO: Add WCS support
        #val = 'X: {:.2f}, Y: {:.2f}'.format(data_x + self._pixel_offset,
        #                                    data_y + self._pixel_offset)
        #val += ', value: {}'.format(imval)
        #self._jup_coord.value = val

    def load_array(self, arr):
        """Load a 2D array into the viewer.

        Parameters
        ----------
        arr : array-like
            2D array.

        """
        x_sc, y_sc = LinearScale(), LinearScale()
        x = np.arange(arr.shape[1])
        y = np.arange(arr.shape[0])
        col_sc = ColorScale(scheme='Greys')
        aspect_ratio = arr.shape[1] / arr.shape[0]
        img = HeatMap(x=x, y=y, color=arr,
                      scales={'x': x_sc, 'y': y_sc, 'color': col_sc})
        ax_x = Axis(scale=x_sc)
        ax_y = Axis(scale=y_sc, orientation='vertical')
        ax_c = ColorAxis(scale=col_sc)

        # FIXME: Does not work
        # https://github.com/bqplot/bqplot/issues/1269
        img.on_hover(self._mouse_move_cb)

        self._jup_img.marks = [img]
        self._jup_img.axes = [ax_x, ax_y, ax_c]
        self._jup_img.min_aspect_ratio = aspect_ratio
        self._jup_img.max_aspect_ratio = aspect_ratio
