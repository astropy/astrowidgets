"""``astrowidgets`` with ``bqplot`` as backend."""

import math
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

# Jupyter widgets
import ipywidgets as ipyw

# bqplot
from bqplot import Figure, HeatMap, LinearScale, ColorScale, Axis, ColorAxis
from bqplot_image_gl.interacts import MouseInteraction

__all__ = ['ImageWidget']

# Allowed locations for cursor display
ALLOWED_CURSOR_LOCATIONS = ['top', 'bottom', None]

# List of marker names that are for internal use only
RESERVED_MARKER_SET_NAMES = ['all']


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
        self._wcs = None

        self._cursor = 'bottom'
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

    def _mouse_cb(self, interaction, data, buffers):
        """
        Callback for mouse events; e.g., to display position in RA/DEC deg.
        """
        event = data['event']
        if event not in ('mousemove', 'click', 'dblclick', 'contextmenu'):
            return  # no-op

        # https://github.com/glue-viz/bqplot-image-gl/pull/37
        if data['event'] == 'mousemove':
            image = self._jup_img.marks[0]
            domain_x = data['domain']['x']
            domain_y = data['domain']['y']
            pixel_x = (domain_x - image.x[0]) / (image.x[1] - image.x[0])
            pixel_y = (domain_y - image.y[0]) / (image.y[1] - image.y[0])
            # TODO: think about +/-1 and pixel edges
            ix = int(math.floor(pixel_x))  # TODO: int(pixel_x + 0.5) ???
            iy = int(math.floor(pixel_y))  # TODO: int(pixel_y + 0.5) ???
            if pixel_x >= 0 and pixel_x < image.color.shape[1] and pixel_y >= 0 and pixel_y < image.color.shape[0]:
                imval = image.color[iy, ix]
            else:
                imval = "NA"
            val = f'X: {pixel_x + self._pixel_offset:.2f}, Y: {pixel_y + self._pixel_offset:.2f}'
            if self._wcs is not None:
                sky = self._wcs.pixel_to_world(pixel_x, pixel_y)
                ra = sky.icrs.ra.to_string(sep='hms')
                dec = sky.icrs.dec.to_string(sep='hms')
                val += f' (RA: {ra}, DEC: {dec})'
            val += f', value: {imval}'

        # TODO: Handle clicks
        else:  # 'click', 'dblclick', 'contextmenu'
            val = f'DEBUG {event}: {data}'

        self._jup_coord.value = val

    def load_fits(self, fitsorfn, numhdu=0, memmap=True):
        """
        Load a FITS file into the viewer.

        Parameters
        ----------
        fitsorfn : str or HDU
            Either a file name or an HDU (*not* an HDUList).
            If file name is given, WCS in primary header is automatically
            inherited. If a single HDU is given, WCS must be in the HDU
            header.

        numhdu : int
            Extension number of the desired HDU.

        memmap : bool
            Memory mapping.

        """
        if isinstance(fitsorfn, str):
            with fits.open(fitsorfn, memmap=memmap) as pf:
                data = pf[numhdu].data
                hdr = pf[numhdu].header
        elif isinstance(fitsorfn, (fits.ImageHDU, fits.CompImageHDU,
                                   fits.PrimaryHDU)):
            data = fitsorfn.data
            hdr = fitsorfn.header
        else:
            return  # no-op

        try:
            wcs = WCS(hdr)
        except Exception as e:
            warnings.warn(f'Failed to load WCS: {repr(e)}')
        else:
            self._wcs = wcs

        self.load_array(data)

    def load_nddata(self, nddata):
        """
        Load an ``NDData`` object into the viewer.

        .. todo:: Add flag/masking support, etc.

        Parameters
        ----------
        nddata : `~astropy.nddata.NDData`
            ``NDData`` with image data and WCS.

        """
        self._wcs = nddata.wcs
        self.load_array(nddata.data)

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

        # TODO: Unset num_ticks after this issue is resolved
        # https://github.com/bqplot/bqplot/issues/1274
        ax_c = ColorAxis(scale=col_sc, num_ticks=5)

        self._jup_img.marks = (img, )
        self._jup_img.axes = [ax_x, ax_y, ax_c]
        self._jup_img.max_aspect_ratio = aspect_ratio
        self._jup_img.min_aspect_ratio = aspect_ratio

        self._jup_img.interaction = MouseInteraction(
            x_scale=img.scales['x'], y_scale=img.scales['y'], move_throttle=70)
        self._jup_img.interaction.on_msg(self._mouse_cb)

    def center_on(self, point):
        """
        Centers the view on a particular point.

        Parameters
        ----------
        point : tuple or `~astropy.coordinates.SkyCoord`
            If tuple of ``(X, Y)`` is given, it is assumed
            to be in data coordinates.
        """
        raise NotImplementedError

    def offset_to(self, dx, dy, skycoord_offset=False):
        """
        Move the center to a point that is given offset
        away from the current center.

        Parameters
        ----------
        dx, dy : float
            Offset value. Unit is assumed based on
            ``skycoord_offset``.

        skycoord_offset : bool
            If `True`, offset must be given in degrees.
            Otherwise, they are in pixel values.

        """
        raise NotImplementedError

    @property
    def zoom_level(self):
        """
        Zoom level:

        * 1 means real-pixel-size.
        * 2 means zoomed in by a factor of 2.
        * 0.5 means zoomed out by a factor of 2.

        """
        raise NotImplementedError

    @zoom_level.setter
    def zoom_level(self, val):
        raise NotImplementedError

    def zoom(self, val):
        """
        Zoom in or out by the given factor.

        Parameters
        ----------
        val : int
            The zoom level to zoom the image.
            See `zoom_level`.

        """
        raise NotImplementedError

    @property
    def is_marking(self):
        """
        `True` if in marking mode, `False` otherwise.
        Marking mode means a mouse click adds a new marker.
        This does not affect :meth:`add_markers`.
        """
        return self._is_marking

    def start_marking(self, marker_name=None,
                      marker=None):
        """
        Start marking, with option to name this set of markers or
        to specify the marker style.
        """
        raise NotImplementedError

    def stop_marking(self, clear_markers=False):
        """
        Stop marking mode, with option to clear markers, if desired.

        Parameters
        ----------
        clear_markers : bool, optional
            If ``clear_markers`` is `False`, existing markers are
            retained until :meth:`reset_markers` is called.
            Otherwise, they are erased.
        """
        raise NotImplementedError

    @property
    def marker(self):
        """
        Marker to use.

        .. todo:: Add more examples.

        Marker can be set as follows::

            {'type': 'circle', 'color': 'cyan', 'radius': 20}
            {'type': 'cross', 'color': 'green', 'radius': 20}
            {'type': 'plus', 'color': 'red', 'radius': 20}

        """
        # Change the marker from a very ginga-specific type (a partial
        # of a ginga drawing canvas type) to a generic dict, which is
        # what we expect the user to provide.
        #
        # That makes things like self.marker = self.marker work.
        return self._marker_dict

    @marker.setter
    def marker(self, val):
        raise NotImplementedError

    def get_markers(self, x_colname='x', y_colname='y',
                    skycoord_colname='coord',
                    marker_name=None):
        """
        Return the locations of existing markers.

        Parameters
        ----------
        x_colname, y_colname : str
            Column names for X and Y data coordinates.
            Coordinates returned are 0- or 1-indexed, depending
            on ``self.pixel_offset``.

        skycoord_colname : str
            Column name for ``SkyCoord``, which contains
            sky coordinates associated with the active image.
            This is ignored if image has no WCS.

        Returns
        -------
        markers_table : `~astropy.table.Table` or ``None``
            Table of markers, if any, or ``None``.

        """
        raise NotImplementedError

    def _validate_marker_name(self, marker_name):
        """
        Raise an error if the marker_name is not allowed.
        """
        raise NotImplementedError

    def add_markers(self, table, x_colname='x', y_colname='y',
                    skycoord_colname='coord', use_skycoord=False,
                    marker_name=None):
        """
        Creates markers in the image at given points.

        .. todo::

            Later enhancements to include more columns
            to control size/style/color of marks,

        Parameters
        ----------
        table : `~astropy.table.Table`
            Table containing marker locations.

        x_colname, y_colname : str
            Column names for X and Y.
            Coordinates can be 0- or 1-indexed, as
            given by ``self.pixel_offset``.

        skycoord_colname : str
            Column name with ``SkyCoord`` objects.

        use_skycoord : bool
            If `True`, use ``skycoord_colname`` to mark.
            Otherwise, use ``x_colname`` and ``y_colname``.

        marker_name : str, optional
            Name to assign the markers in the table. Providing a name
            allows markers to be removed by name at a later time.
        """
        raise NotImplementedError

    def remove_markers(self, marker_name=None):
        """
        Remove some but not all of the markers by name used when
        adding the markers

        Parameters
        ----------

        marker_name : str, optional
            Name used when the markers were added.
        """
        raise NotImplementedError

    def reset_markers(self):
        """
        Delete all markers.
        """
        raise NotImplementedError

    @property
    def stretch_options(self):
        """
        List all available options for image stretching.
        """
        raise NotImplementedError

    @property
    def stretch(self):
        """
        The image stretching algorithm in use.
        """
        raise NotImplementedError

    @stretch.setter
    def stretch(self, val):
        raise NotImplementedError

    @property
    def autocut_options(self):
        """
        List all available options for image auto-cut.
        """
        raise NotImplementedError

    @property
    def cuts(self):
        """
        Current image cut levels.
        To set new cut levels, either provide a tuple of
        ``(low, high)`` values or one of the options from
        `autocut_options`.
        """
        raise NotImplementedError

    # TODO: Possible to use astropy.visualization directly?
    @cuts.setter
    def cuts(self, val):
        raise NotImplementedError

    @property
    def colormap_options(self):
        """List of colormap names."""
        raise NotImplementedError

    def set_colormap(self, cmap):
        """
        Set colormap to the given colormap name.

        Parameters
        ----------
        cmap : str
            Colormap name. Possible values can be obtained from
            :meth:`colormap_options`.

        """
        raise NotImplementedError

    @property
    def cursor(self):
        """
        Show or hide cursor information (X, Y, WCS).
        Acceptable values are 'top', 'bottom', or ``None``.
        """
        return self._cursor

    @cursor.setter
    def cursor(self, val):
        if val is None:
            self._jup_coord.layout.visibility = 'hidden'
            self._jup_coord.layout.display = 'none'
        elif val == 'top' or val == 'bottom':
            self._jup_coord.layout.visibility = 'visible'
            self._jup_coord.layout.display = 'flex'
            if val == 'top':
                self.layout.flex_flow = 'column-reverse'
            else:
                self.layout.flex_flow = 'column'
        else:
            raise ValueError(f'Invalid value {val} for cursor.'
                             f'Valid values are: {ALLOWED_CURSOR_LOCATIONS}')
        self._cursor = val

    @property
    def click_center(self):
        """
        Settable.
        If True, middle-clicking can be used to center.  If False, that
        interaction is disabled.

        In the future this might go from True/False to being a selectable
        button. But not for the first round.
        """
        raise NotImplementedError

    @click_center.setter
    def click_center(self, val):
        raise NotImplementedError

    @property
    def click_drag(self):
        """
        Settable.
        If True, the "click-and-drag" mode is an available interaction for
        panning.  If False, it is not.

        Note that this should be automatically made `False` when selection mode
        is activated.
        """
        raise NotImplementedError

    @click_drag.setter
    def click_drag(self, value):
        raise NotImplementedError

    @property
    def scroll_pan(self):
        """
        Settable.
        If True, scrolling moves around in the image.  If False, scrolling
        (up/down) *zooms* the image in and out.
        """
        raise NotImplementedError

    @scroll_pan.setter
    def scroll_pan(self, value):
        raise NotImplementedError

    def save(self, filename):
        """
        Save out the current image view to given PNG filename.
        """
        raise NotImplementedError
