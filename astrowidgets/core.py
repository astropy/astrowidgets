"""Module containing core functionality of ``astrowidgets``."""

# STDLIB
import functools
import warnings

# THIRD-PARTY
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table

# Jupyter widgets
import ipywidgets as ipyw

# Ginga
from ginga.AstroImage import AstroImage
from ginga.canvas.CanvasObject import drawCatalog
from ginga.web.jupyterw.ImageViewJpw import EnhancedCanvasView
from ginga.util.wcs import raDegToString, decDegToString

__all__ = ['ImageWidget']

# Allowed locations for cursor display
ALLOWED_CURSOR_LOCATIONS = ['top', 'bottom', None]


class ImageWidget(ipyw.VBox):
    """
    Image widget for Jupyter notebook using Ginga viewer.

    .. todo:: Any property passed to constructor has to be valid keyword.

    Parameters
    ----------
    logger : obj or ``None``
        Ginga logger. For example::

            from ginga.misc.log import get_logger
            logger = get_logger('my_viewer', log_stderr=False,
                                log_file='ginga.log', level=40)

    image_width, image_height : int
        Dimension of Jupyter notebook's image widget.

    use_opencv : bool
        Let Ginga use ``opencv`` to speed up image transformation;
        e.g., rotation and mosaic. If this is enabled and you
        do not have ``opencv``, you will get a warning.

    """

    def __init__(self, logger=None, image_width=500, image_height=500,
                 use_opencv=True):
        super().__init__()

        # TODO: Is this the best place for this?
        if use_opencv:
            try:
                from ginga import trcalc
                trcalc.use('opencv')
            except ImportError as exc:
                warnings.warn('install opencv or set use_opencv=False')

        self._viewer = EnhancedCanvasView(logger=logger)
        self._is_marking = False
        self._click_center = False

        self._jup_img = ipyw.Image(format='jpeg')

        # Set the image margin to over the widgets default of 2px on
        # all sides.
        self._jup_img.layout.margin = '0'

        # Set both of those to ensure consistent display in notebook
        # and jupyterlab when the image is put into a container smaller
        # than the image.

        self._jup_img.max_width = '100%'
        self._jup_img.height = 'auto'

        # Set the width of the box containing the image to the desired width
        self.layout.width = str(image_width)

        # Note we are NOT setting the height. That is because the height
        # is automatically set by the image aspect ratio.

        # These need to also be set for now; ginga uses them to figure
        # out what size image to make.
        self._jup_img.width = image_width
        self._jup_img.height = image_height

        self._viewer.set_widget(self._jup_img)

        # enable all possible keyboard and pointer operations
        self._viewer.get_bindings().enable_all(True)

        # enable draw
        self.dc = drawCatalog
        self.canvas = self.dc.DrawingCanvas()
        self.canvas.enable_draw(True)
        self.canvas.enable_edit(True)

        # Marker
        self.marker = {'type': 'circle', 'color': 'cyan', 'radius': 20}
        self._marktag = 'marktag'

        # coordinates display
        self._jup_coord = ipyw.HTML('Coordinates show up here')
        # This needs ipyevents 0.3.1 to work
        self._viewer.add_callback('cursor-changed', self._mouse_move_cb)
        self._viewer.add_callback('cursor-down', self._mouse_click_cb)

        # Define a callback that shows the output of a print
        self.print_out = ipyw.Output()

        self._cursor = 'bottom'
        self.children = [self._jup_img, self._jup_coord]

    @property
    def logger(self):
        """Logger for this widget."""
        return self._viewer.logger

    @property
    def image_width(self):
        return int(self._jup_img.width)

    @image_width.setter
    def image_width(self, value):
        # widgets expect width/height as strings, but most users will not, so
        # do the conversion.
        self._jup_img.width = str(value)
        self._viewer.set_window_size(self.image_width, self.image_height)

    @property
    def image_height(self):
        return int(self._jup_img.height)

    @image_height.setter
    def image_height(self, value):
        # widgets expect width/height as strings, but most users will not, so
        # do the conversion.
        self._jup_img.height = str(value)
        self._viewer.set_window_size(self.image_width, self.image_height)

    def _mouse_move_cb(self, viewer, button, data_x, data_y):
        """
        Callback to display position in RA/DEC deg.
        """
        if self.cursor is None:  # no-op
            return

        image = viewer.get_image()
        if image is not None:
            ix = int(data_x + 0.5)
            iy = int(data_y + 0.5)
            try:
                imval = viewer.get_data(ix, iy)
            except Exception:
                imval = 'N/A'

            # Same as setting pixel_coords_offset=1 in general.cfg
            val = 'X: {:.2f}, Y:{:.2f}'.format(data_x + 1, data_y + 1)
            if image.wcs.wcs is not None:
                ra, dec = image.pixtoradec(data_x, data_y)
                val += ' (RA: {}, DEC: {})'.format(
                    raDegToString(ra), decDegToString(dec))
            val += ', value: {}'.format(imval)
            self._jup_coord.value = val

    def _mouse_click_cb(self, viewer, event, data_x, data_y):
        """
        Callback to handle mouse clicks.
        """
        if self.is_marking:
            objs = []
            try:
                c_mark = viewer.canvas.get_object_by_tag(self._marktag)
            except Exception:  # Nothing drawn yet
                pass
            else:  # Add to existing marks
                objs = c_mark.objects
                viewer.canvas.delete_object_by_tag(self._marktag)

            # NOTE: By always using CompoundObject, marker handling logic
            # is simplified.
            obj = self.marker(x=data_x, y=data_y)
            objs.append(obj)
            self._marktag = viewer.canvas.add(self.dc.CompoundObject(*objs))

            with self.print_out:
                print('Selected {} {}'.format(obj.x, obj.y))

        elif self.click_center:
            self.center_on((data_x, data_y))

            with self.print_out:
                print('Centered on X={} Y={}'.format(data_x + 1, data_y + 1))

#     def _repr_html_(self):
#         """
#         Show widget in Jupyter notebook.
#         """
#         from IPython.display import display
#         return display(self._widget)

    def load_fits(self, fitsorfn, numhdu=None, memmap=None):
        """
        Load a FITS file into the viewer.

        Parameters
        ----------
        fitsorfn : str or HDU
            Either a file name or an HDU (*not* an HDUList).
            If file name is given, WCS in primary header is automatically
            inherited. If a single HDU is given, WCS must be in the HDU
            header.

        numhdu : int or ``None``
            Extension number of the desired HDU.
            If ``None``, it is determined automatically.

        memmap : bool or ``None``
            Memory mapping.
            If ``None``, it is determined automatically.

        """
        if isinstance(fitsorfn, str):
            image = AstroImage(logger=self.logger, inherit_primary_header=True)
            image.load_file(fitsorfn, numhdu=numhdu, memmap=memmap)
            self._viewer.set_image(image)

        elif isinstance(fitsorfn, (fits.ImageHDU, fits.CompImageHDU,
                                   fits.PrimaryHDU)):
            self._viewer.load_hdu(fitsorfn)

    def load_nddata(self, nddata):
        """
        Load an ``NDData`` object into the viewer.

        .. todo:: Add flag/masking support, etc.

        Parameters
        ----------
        nddata : `~astropy.nddata.NDData`
            ``NDData`` with image data and WCS.

        """
        from ginga.util.wcsmod.wcs_astropy import AstropyWCS

        image = AstroImage(logger=self.logger)
        image.set_data(nddata.data)
        _wcs = AstropyWCS(self.logger)
        if nddata.wcs:
            _wcs.load_header(nddata.wcs.to_header())

        try:
            image.set_wcs(_wcs)
        except Exception as e:
            print('Unable to set WCS from NDData: {}'.format(str(e)))
        self._viewer.set_image(image)

    def load_array(self, arr):
        """
        Load a 2D array into the viewer.

        .. note:: Use :meth:`load_nddata` for WCS support.

        Parameters
        ----------
        arr : array-like
            2D array.

        """
        self._viewer.load_data(arr)

    def center_on(self, point, pixel_coords_offset=1):
        """
        Centers the view on a particular point.

        Parameters
        ----------
        point : tuple or `~astropy.coordinates.SkyCoord`
            If tuple of ``(X, Y)`` is given, it is assumed
            to be in data coordinates.

        pixel_coords_offset : {0, 1}
            Data coordinates provided are n-indexed,
            where n is the given value.
            This is ignored if ``SkyCoord`` is provided.

        """
        if isinstance(point, SkyCoord):
            self._viewer.set_pan(point.ra.deg, point.dec.deg, coord='wcs')
        else:
            self._viewer.set_pan(*(np.asarray(point) - pixel_coords_offset))

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
        if skycoord_offset:
            coord = 'wcs'
        else:
            coord = 'data'

        pan_x, pan_y = self._viewer.get_pan(coord=coord)
        self._viewer.set_pan(pan_x + dx, pan_y + dy, coord=coord)

    @property
    def zoom_level(self):
        """
        Zoom level:

        * 1 means real-pixel-size.
        * 2 means zoomed in by a factor of 2.
        * 0.5 means zoomed out by a factor of 2.

        """
        return self._viewer.get_scale()

    @zoom_level.setter
    def zoom_level(self, val):
        self._viewer.scale_to(val, val)

    def zoom(self, val):
        """
        Zoom in or out by the given factor.

        Parameters
        ----------
        val : int
            The zoom level to zoom the image.
            See `zoom_level`.

        """
        self.zoom_level = self.zoom_level * val

    @property
    def is_marking(self):
        """
        `True` if in marking mode, `False` otherwise.
        Marking mode means a mouse click adds a new marker.
        This does not affect :meth:`add_markers`.
        """
        return self._is_marking

    @is_marking.setter
    def is_marking(self, val):
        if not isinstance(val, bool):
            raise ValueError('Must be True or False')
        elif self.click_center and val:
            raise ValueError('Cannot set to True while in click-center mode')
        self._is_marking = val

    def stop_marking(self, clear_markers=True):
        """
        Stop marking mode, with option to clear markers, if desired.

        Parameters
        ----------
        clear_markers : bool
            If ``clear_markers`` is `False`, existing markers are
            retained until :meth:`reset_markers` is called.
            Otherwise, they are erased.
        """
        self.is_marking = False
        if clear_markers:
            self.reset_markers()

    @property
    def marker(self):
        """
        Marker to use.

        .. todo:: Add more examples.

        Marker can be set as follows::

            {'type': 'circle', 'color': 'cyan', 'radius': 20}

        """
        return self._marker

    @marker.setter
    def marker(self, val):
        marker_type = val.pop('type')
        if marker_type == 'circle':
            self._marker = functools.partial(self.dc.Circle, **val)
        else:  # TODO: Implement more shapes
            raise NotImplementedError(
                'Marker type "{}" not supported'.format(marker_type))

    def get_markers(self, x_colname='x', y_colname='y',
                    pixel_coords_offset=1,
                    skycoord_colname='coord'):
        """
        Return the locations of existing markers.

        Parameters
        ----------
        x_colname, y_colname : str
            Column names for X and Y data coordinates.
            Coordinates returned are 0- or 1-indexed, depending
            on ``pixel_coords_offset``.

        pixel_coords_offset : {0, 1}
            Data coordinates returned are n-indexed,
            where n is the given value.

        skycoord_colname : str
            Column name for ``SkyCoord``, which contains
            sky coordinates associated with the active image.
            This is ignored if image has no WCS.

        Returns
        -------
        markers_table : `~astropy.table.Table` or ``None``
            Table of markers, if any, or ``None``.

        """
        try:
            c_mark = self._viewer.canvas.get_object_by_tag(self._marktag)
        except Exception as e:  # No markers
            self.logger.warning(str(e))
            return

        image = self._viewer.get_image()
        xy_col = []

        if image.wcs.wcs is None:  # Do not include SkyCoord column
            include_skycoord = False
        else:
            include_skycoord = True
            radec_col = []

        # Extract coordinates from markers
        for obj in c_mark.objects:
            if obj.coord == 'data':
                xy_col.append([obj.x, obj.y])
                if include_skycoord:
                    radec_col.append([np.nan, np.nan])
            elif not include_skycoord:  # marker in WCS but image has none
                self.logger.warning(
                    'Skipping ({},{}); image has no WCS'.format(obj.x, obj.y))
            else:  # wcs
                xy_col.append([np.nan, np.nan])
                radec_col.append([obj.x, obj.y])

        # Convert to numpy arrays
        xy_col = np.asarray(xy_col)  # [[x0, y0], [x1, y1], ...]

        if include_skycoord:
            # [[ra0, dec0], [ra1, dec1], ...]
            radec_col = np.asarray(radec_col)

            # Fill in X,Y from RA,DEC
            mask = np.isnan(xy_col[:, 0])  # One bool per row
            if np.any(mask):
                xy_col[mask] = image.wcs.wcspt_to_datapt(radec_col[mask])

            # Fill in RA,DEC from X,Y
            mask = np.isnan(radec_col[:, 0])
            if np.any(mask):
                radec_col[mask] = image.wcs.datapt_to_wcspt(xy_col[mask])

            sky_col = SkyCoord(radec_col[:, 0], radec_col[:, 1], unit='deg')

        # Convert X,Y from 0-indexed to 1-indexed
        if pixel_coords_offset != 0:
            xy_col += pixel_coords_offset

        # Build table
        if include_skycoord:
            markers_table = Table(
                [xy_col[:, 0], xy_col[:, 1], sky_col],
                names=(x_colname, y_colname, skycoord_colname))
        else:
            markers_table = Table(xy_col.T, names=(x_colname, y_colname))

        return markers_table

    def add_markers(self, table, x_colname='x', y_colname='y',
                    pixel_coords_offset=1,
                    skycoord_colname='coord', use_skycoord=False):
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
            given by ``pixel_coords_offset``.

        pixel_coords_offset : {0, 1}
            Data coordinates provided are n-indexed,
            where n is the given value.
            This is ignored if ``use_skycoord=True``.

        skycoord_colname : str
            Column name with ``SkyCoord`` objects.

        use_skycoord : bool
            If `True`, use ``skycoord_colname`` to mark.
            Otherwise, use ``x_colname`` and ``y_colname``.

        """
        # TODO: Resolve https://github.com/ejeschke/ginga/issues/672
        # Extract coordinates from table.
        # They are always arrays, not scalar.
        if use_skycoord:
            image = self._viewer.get_image()
            if image is None:
                raise ValueError('Cannot get image from viewer')
            if image.wcs.wcs is None:
                raise ValueError(
                    'Image has no valid WCS, '
                    'try again with use_skycoord=False')
            coord_type = 'wcs'
            coord_val = table[skycoord_colname]
            coord_x = coord_val.ra.deg
            coord_y = coord_val.dec.deg
        else:  # Use X,Y
            coord_type = 'data'
            coord_x = table[x_colname].data
            coord_y = table[y_colname].data
            # Convert data coordinates from 1-indexed to 0-indexed
            if pixel_coords_offset != 0:
                coord_x -= pixel_coords_offset
                coord_y -= pixel_coords_offset

        # Prepare canvas and retain existing marks
        objs = []
        try:
            c_mark = self._viewer.canvas.get_object_by_tag(self._marktag)
        except Exception:
            pass
        else:
            objs = c_mark.objects
            self._viewer.canvas.delete_object_by_tag(self._marktag)

        # TODO: Test to see if we can mix WCS and data on the same canvas
        objs += [self.marker(x=x, y=y, coord=coord_type)
                 for x, y in zip(coord_x, coord_y)]
        self._marktag = self._viewer.canvas.add(self.dc.CompoundObject(*objs))

    # TODO: Future work?
    def remove_markers(self, arr):
        """
        Remove some but not all of the markers.

        Parameters
        ----------
        arr : ``SkyCoord`` or array-like
            Sky coordinates or 2xN array.

        """
        # NOTE: How to match? Use np.isclose?
        #       What if there are 1-to-many matches?
        return NotImplementedError

    def reset_markers(self):
        """
        Delete all markers.
        """
        try:
            self._viewer.canvas.delete_object_by_tag(self._marktag)
        except Exception:
            pass

    @property
    def stretch_options(self):
        """
        List all available options for image stretching.
        """
        return self._viewer.get_color_algorithms()

    @property
    def stretch(self):
        """
        The image stretching algorithm in use.
        """
        return self._viewer.rgbmap.dist

    # TODO: Possible to use astropy.visualization directly?
    @stretch.setter
    def stretch(self, val):
        valid_vals = self.stretch_options
        if val not in valid_vals:
            raise ValueError('Value must be one of: {}'.format(valid_vals))
        self._viewer.set_color_algorithm(val)

    @property
    def autocut_options(self):
        """
        List all available options for image auto-cut.
        """
        return self._viewer.get_autocut_methods()

    @property
    def cuts(self):
        """
        Current image cut levels.
        To set new cut levels, either provide a tuple of
        ``(low, high)`` values or one of the options from
        `autocut_options`.
        """
        return self._viewer.get_cut_levels()

    # TODO: Possible to use astropy.visualization directly?
    @cuts.setter
    def cuts(self, val):
        if isinstance(val, str):  # Autocut
            valid_vals = self.autocut_options
            if val not in valid_vals:
                raise ValueError('Value must be one of: {}'.format(valid_vals))
            self._viewer.set_autocut_params(val)
        else:  # (low, high)
            if len(val) > 2:
                raise ValueError('Value must have length 2.')
            self._viewer.cut_levels(val[0], val[1])

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
            raise ValueError('Invalid value {} for cursor.'.format(val) +
                             'Valid values are: '
                             '{}'.format(ALLOWED_CURSOR_LOCATIONS))
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
        return self._click_center

    @click_center.setter
    def click_center(self, val):
        if not isinstance(val, bool):
            raise ValueError('Must be True or False')
        elif self.is_marking and val:
            raise ValueError('Cannot set to True while in marking mode')
        self._click_center = val

    # TODO: Awaiting https://github.com/ejeschke/ginga/issues/674
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

    @property
    def scroll_pan(self):
        """
        Settable.
        If True, scrolling moves around in the image.  If False, scrolling
        (up/down) *zooms* the image in and out.
        """
        raise NotImplementedError

    # https://github.com/ejeschke/ginga/pull/665
    def save(self, filename):
        """
        Save out the current image view to given PNG filename.
        """
        self._viewer.save_rgb_image_as_file(filename)
