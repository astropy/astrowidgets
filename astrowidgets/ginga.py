"""The ``astrowidgets.ginga`` module contains a widget implemented with the
Ginga backend.

For this to work, ``astrowidgets`` must be installed along with the optional
dependencies specified for the Ginga backend; e.g.,::

    pip install 'astrowidgets[ginga]'

"""
import functools
import os
import warnings

import ipywidgets as ipyw
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.utils.decorators import deprecated

from ginga.AstroImage import AstroImage
from ginga.canvas.CanvasObject import drawCatalog
from ginga.web.jupyterw.ImageViewJpw import EnhancedCanvasView
from ginga.util.wcs import ra_deg_to_str, dec_deg_to_str

__all__ = ['ImageWidget']


class ImageWidget(ipyw.VBox):
    """Image widget for Jupyter notebook using Ginga viewer.

    .. todo:: Any property passed to constructor has to be valid keyword.

    Parameters
    ----------
    logger : obj or None
        Ginga logger. For example::

            from ginga.misc.log import get_logger
            logger = get_logger('my_viewer', log_stderr=False,
                                log_file='ginga.log', level=40)

    image_width, image_height : int
        Dimension of Jupyter notebook's image widget.

    image_widget : obj or None
        Jupyter notebook's image widget. If None, a new widget will be created.

    cursor_widget : obj or None
        Jupyter notebook's cursor widget. If None, a new widget will be created.
    """

    def __init__(self, logger=None, image_width=500, image_height=500,
                 image_widget=None, cursor_widget=None,
                 **kwargs):
        super().__init__()
        if 'use_opencv' in kwargs:
            warnings.warn("use_opencv kwarg has been deprecated--"
                          "opencv will be used if it is installed",
                          DeprecationWarning)

        self._viewer = EnhancedCanvasView(logger=logger)
        self.ALLOWED_CURSOR_LOCATIONS = ['top', 'bottom', None]
        self.RESERVED_MARKER_SET_NAMES = ['all']

        if image_widget is None:
            self._jup_img = ipyw.Image(format='jpeg')
        else:
            self._jup_img = image_widget

        if cursor_widget is None:
            self._jup_coord = ipyw.HTML('Coordinates show up here')
        else:
            self._jup_coord = cursor_widget

        if isinstance(self._jup_img, ipyw.Image):
            # Set the image margin on all sides.
            self._jup_img.layout.margin = '0'

            # Set both of those to ensure consistent display in notebook
            # and jupyterlab when the image is put into a container smaller
            # than the image.
            self._jup_img.max_width = '100%'
            self._jup_img.height = 'auto'

        # Set the width of the box containing the image to the desired width
        # Note: We are NOT setting the height. That is because the height
        # is automatically set by the image aspect ratio.
        self.layout.width = str(image_width)

        # Make sure all of the internal state trackers have a value
        # and start in a state which is definitely allowed: all are
        # False.
        self._is_marking = False
        self._click_center = False
        self._click_drag = False
        self._scroll_pan = False
        self._cached_state = {}

        # Marker
        self._marker_dict = {}
        self._marker = None
        # Maintain marker tags as a set because we do not want
        # duplicate names.
        self._marktags = set()
        # Let's have a default name for the tag too:
        self._default_mark_tag_name = 'default-marker-name'
        self._interactive_marker_set_name_default = 'interactive-markers'
        self._interactive_marker_set_name = self._interactive_marker_set_name_default

        # Define a callback that shows the output of a print
        self.print_out = ipyw.Output()

        self._cursor = 'bottom'
        self.children = [self._jup_img, self._jup_coord]

        # These need to also be set for now.
        # Ginga uses them to figure out what size image to make.
        self._jup_img.width = image_width
        self._jup_img.height = image_height

        self._viewer = EnhancedCanvasView(logger=logger)
        self._viewer.set_widget(self._jup_img)

        # Enable all possible keyboard and pointer operations
        self._viewer.get_bindings().enable_all(True)

        # Enable draw
        self.dc = drawCatalog
        self.canvas = self.dc.DrawingCanvas()
        self.canvas.enable_draw(True)
        self.canvas.enable_edit(True)

        # Set a couple of things to match the Ginga defaults
        self._scroll_pan = True
        self._click_drag = False

        bind_map = self._viewer.get_bindmap()
        # Set up right-click and drag adjusts the contrast
        bind_map.map_event(None, (), 'ms_right', 'contrast')
        # Shift-right-click restores the default contrast
        bind_map.map_event(None, ('shift',), 'ms_right', 'contrast_restore')

        # Marker
        self.marker = {'type': 'circle', 'color': 'cyan', 'radius': 20}

        # This needs ipyevents 0.3.1 to work
        self._viewer.add_callback('cursor-changed', self._mouse_move_cb)
        self._viewer.add_callback('cursor-down', self._mouse_click_cb)

    @property
    def viewer(self):
        return self._viewer

    @property
    def logger(self):
        """Logger for this widget."""
        return self._viewer.logger

    # Need this here because we need to overwrite the setter.
    @property
    def image_width(self):
        """Width of image widget."""
        return int(self._jup_img.width)

    @image_width.setter
    def image_width(self, value):
        # widgets expect width/height as strings, but most users will not, so
        # do the conversion.
        self._jup_img.width = str(value)
        self._viewer.set_window_size(self.image_width, self.image_height)

    # Need this here because we need to overwrite the setter.
    @property
    def image_height(self):
        """Height of image widget."""
        return int(self._jup_img.height)

    @image_height.setter
    def image_height(self, value):
        # widgets expect width/height as strings, but most users will not, so
        # do the conversion.
        self._jup_img.height = str(value)
        self._viewer.set_window_size(self.image_width, self.image_height)

    def _mouse_move_cb(self, viewer, button, data_x, data_y):
        """Callback to display position in RA/DEC deg."""
        if self.cursor is None:  # no-op
            return

        image = viewer.get_image()
        if image is not None:
            ix = int(data_x + 0.5)
            iy = int(data_y + 0.5)
            try:
                imval = viewer.get_data(ix, iy)
                imval = f'{imval:8.3f}'
            except Exception:
                imval = 'N/A'

            val = (f'X: {data_x:.2f}, '
                   f'Y: {data_y:.2f}')

            if image.wcs.wcs is not None:
                try:
                    ra, dec = image.pixtoradec(data_x, data_y)
                    val += (f' (RA: {ra_deg_to_str(ra)},'
                            f' DEC: {dec_deg_to_str(dec)})')
                except Exception:
                    val += ' (RA, DEC: WCS error)'

            val += f', value: {imval}'
            self._jup_coord.value = val

    def _mouse_click_cb(self, viewer, event, data_x, data_y):
        """Callback to handle mouse clicks."""
        if self.is_marking:
            marker_name = self._interactive_marker_set_name
            objs = []
            try:
                c_mark = viewer.canvas.get_object_by_tag(marker_name)
            except Exception:  # Nothing drawn yet
                pass
            else:  # Add to existing marks
                objs = c_mark.objects
                viewer.canvas.delete_object_by_tag(marker_name)

            # NOTE: By always using CompoundObject, marker handling logic
            # is simplified.
            obj = self._marker(x=data_x, y=data_y, coord='data')
            objs.append(obj)
            viewer.canvas.add(self.dc.CompoundObject(*objs), tag=marker_name)
            self._marktags.add(marker_name)

            # For debugging.
            with self.print_out:
                print(f'Selected {obj.x} {obj.y}')

        elif self.click_center:
            self.center_on((data_x, data_y))

            # For debugging.
            with self.print_out:
                print(f'Centered on X={data_x} '
                      f'Y={data_y}')

    def load_fits(self, filename, numhdu=None, memmap=None):
        """Load a FITS file or HDU into the viewer.

        Parameters
        ----------
        filename : str or HDU
            Name of the FITS file or a HDU (*not* a ``HDUList``).
            If a filename is given, any information in the primary header,
            including WCS, is automatically inherited. If a HDU is given,
            the WCS must be in the HDU header.

        numhdu : int or `None`
            Extension number of the desired HDU. If not given, it is
            determined automatically. This is only used if a filename is given.

        memmap : bool or `None`
            Memory mapping. See `astropy.io.fits.open`.
            This is only used if a filename is given.

        Raises
        ------
        ValueError
           Given ``filename`` type is not supported.

        """
        if isinstance(filename, str):
            image = AstroImage(logger=self.logger, inherit_primary_header=True)
            image.load_file(filename, numhdu=numhdu, memmap=memmap)
            self._viewer.set_image(image)

        elif isinstance(filename, (fits.ImageHDU, fits.CompImageHDU,
                                   fits.PrimaryHDU)):
            self._viewer.load_hdu(filename)
        else:
            raise ValueError(f'Unable to open {filename}')

    def load_nddata(self, nddata):
        """Load a `~astropy.nddata.NDData` object into the viewer.

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
            self.logger.warning(f'Unable to set WCS from NDData: {repr(e)}')
        self._viewer.set_image(image)

    def load_array(self, arr):
        """Load a 2D array into the viewer.

        .. note:: Use :meth:`load_nddata` for WCS support.

        Parameters
        ----------
        arr : array-like
            2D array.

        """
        self._viewer.load_data(arr)

    def center_on(self, point):
        if isinstance(point, SkyCoord):
            self._viewer.set_pan(point.ra.deg, point.dec.deg, coord='wcs')
        else:
            self._viewer.set_pan(*(np.asarray(point)))

    @deprecated('0.3', alternative='offset_by')
    def offset_to(self, dx, dy, skycoord_offset=False):
        """
        Move the center to a point that is given offset
        away from the current center.

        .. note:: This is deprecated. Use :meth:`offset_by`.

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

    def offset_by(self, dx, dy):
        """
        Move the center to a point that is given offset
        away from the current center.

        Parameters
        ----------
        dx, dy : float or `~astropy.unit.Quantity`
            Offset value. Without a unit, assumed to be pixel offsets.
            If a unit is attached, offset by pixel or sky is assumed from
            the unit.

        """
        dx_val, dx_coord = _offset_is_pixel_or_sky(dx)
        dy_val, dy_coord = _offset_is_pixel_or_sky(dy)

        if dx_coord != dy_coord:
            raise ValueError(f'dx is of type {dx_coord} but dy is of type {dy_coord}')

        pan_x, pan_y = self._viewer.get_pan(coord=dx_coord)
        self._viewer.set_pan(pan_x + dx_val, pan_y + dy_val, coord=dx_coord)

    @property
    def zoom_level(self):
        """Zoom level (settable):

        * 1 means real-pixel-size.
        * 2 means zoomed in by a factor of 2.
        * 0.5 means zoomed out by a factor of 2.
        * ``fit`` means fit the image to the viewer.
        """
        return self._viewer.get_scale()

    @zoom_level.setter
    def zoom_level(self, value):
        if value == 'fit':
            self._viewer.zoom_fit()
        else:
            self._viewer.scale_to(value, value)


    def zoom(self, value):
        """Zoom in or out by the given factor.

        Parameters
        ----------
        value : int
            The zoom level to zoom the image.
            See `zoom_level`.

        """
        self.zoom_level = self.zoom_level * value

    @property
    def is_marking(self):
        """`True` if in marking mode, `False` otherwise.
        Marking mode means a mouse click adds a new marker.
        This does not affect :meth:`add_markers`.

        """
        return self._is_marking

    def start_marking(self, marker_name=None, marker=None):
        """Start marking, with option to name this set of markers or
        to specify the marker style.

        This disables `click_center` and `click_drag`, but enables `scroll_pan`.

        Parameters
        ----------
        marker_name : str or `None`, optional
            Marker name to use. This is useful if you want to set different
            groups of markers. If given, this cannot be already defined in
            ``RESERVED_MARKER_SET_NAMES`` attribute. If not given, an internal
            default is used.

        marker : dict or `None`, optional
            Set the marker properties; see `marker`. If not given, the current
            setting is used.

        """
        self.set_cached_state()
        self.click_center = False
        self.click_drag = False
        self.scroll_pan = True  # Set this to ensure there is a mouse way to pan
        self._is_marking = True
        if marker_name is not None:
            self.validate_marker_name(marker_name)
            self._interactive_marker_set_name = marker_name
            self._marktags.add(marker_name)
        else:
            self._interactive_marker_set_name = self._interactive_marker_set_name_default
        if marker is not None:
            self.marker = marker

    def stop_marking(self, clear_markers=False):
        """Stop marking mode, with option to clear all markers, if desired.

        Parameters
        ----------
        clear_markers : bool, optional
            If `False`, existing markers are retained until
            :meth:`remove_all_markers` is called.
            Otherwise, they are all erased.

        """
        if self.is_marking:
            self._is_marking = False
            self.restore_and_clear_cached_state()
            if clear_markers:
                self.remove_all_markers()

    @property
    def marker(self):
        """A dictionary defining the current marker properties.

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
    def marker(self, value):
        # Make a new copy to avoid modifying the dict that the user passed in.
        _marker = value.copy()
        marker_type = _marker.pop('type')
        if marker_type == 'circle':
            self._marker = functools.partial(self.dc.Circle, **_marker)
        elif marker_type == 'plus':
            _marker['type'] = 'point'
            _marker['style'] = 'plus'
            self._marker = functools.partial(self.dc.Point, **_marker)
        elif marker_type == 'cross':
            _marker['type'] = 'point'
            _marker['style'] = 'cross'
            self._marker = functools.partial(self.dc.Point, **_marker)
        else:  # TODO: Implement more shapes
            raise NotImplementedError(f'Marker type "{marker_type}" not supported')
        # Only set this once we have successfully created a marker
        self._marker_dict = value

    def get_marker_names(self):
        """Return a list of used marker names.

        Returns
        -------
        names : list of str
            Sorted list of marker names.

        """
        return sorted(self._marktags)

    def get_markers_by_name(self, marker_name, x_colname='x', y_colname='y',
                            skycoord_colname='coord'):

        # We should always allow the default name. The case
        # where that table is empty will be handled in a moment.
        if (marker_name not in self._marktags
                and marker_name != self._default_mark_tag_name):
            raise ValueError(f"No markers named '{marker_name}' found.")

        try:
            c_mark = self._viewer.canvas.get_object_by_tag(marker_name)
        except Exception:
            # No markers in this table. Issue a warning and continue.
            # Test wants this outside of logger, so...
            warnings.warn(f"Marker set named '{marker_name}' is empty", UserWarning)
            return None

        image = self._viewer.get_image()
        xy_col = []

        if (image is None) or (image.wcs.wcs is None):
            # Do not include SkyCoord column
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
            elif not include_skycoord:  # Marker in WCS but image has none
                self.logger.warning(f'Skipping ({obj.x},{obj.y}); image has no WCS')
            else:  # WCS
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

        # Build table
        if include_skycoord:
            markers_table = Table(
                [xy_col[:, 0], xy_col[:, 1], sky_col],
                names=(x_colname, y_colname, skycoord_colname))
        else:
            markers_table = Table(xy_col, names=(x_colname, y_colname))

        # Either way, add the marker names
        markers_table['marker name'] = marker_name
        return markers_table


    def get_all_markers(self, x_colname='x', y_colname='y', skycoord_colname='coord'):
        """Run :meth:`get_markers_by_name` for all markers."""

        # If it wasn't for the fact that SkyCoord columns can't
        # be stacked this would all fit nicely into a list
        # comprehension. But they can't, so we delete the
        # SkyCoord column if it is present, then add it
        # back after we have stacked.
        coordinates = []
        tables = []
        for name in self._marktags:
            table = self.get_markers_by_name(
                name, x_colname=x_colname, y_colname=y_colname,
                skycoord_colname=skycoord_colname)
            if table is None:
                continue  # No markers by this name, skip it

            if skycoord_colname in table.colnames:
                coordinates.extend(c for c in table[skycoord_colname])
                del table[skycoord_colname]

            tables.append(table)

        if len(tables) == 0:
            return None

        stacked = vstack(tables, join_type='exact')

        if coordinates:
            n_rows = len(stacked)
            n_coo = len(coordinates)
            if n_coo != n_rows:  # This guards against Table auto-broadcast
                raise ValueError(f'Expects {n_rows} coordinates but found {n_coo},'
                                 'some markers may be corrupted')
            stacked[skycoord_colname] = SkyCoord(coordinates)

        return stacked

    # TODO: Resolve https://github.com/ejeschke/ginga/issues/672
    # TODO: Later enhancements to include more columns to control
    # size/style/color of marks
    def add_markers(self, table, x_colname='x', y_colname='y',
                    skycoord_colname='coord', use_skycoord=False,
                    marker_name=None):

        # For now we always convert marker locations to pixels; see
        # comment below.
        coord_type = 'data'

        if marker_name is None:
            marker_name = self._default_mark_tag_name

        self.validate_marker_name(marker_name)
        self._marktags.add(marker_name)

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
            coord_val = table[skycoord_colname]
            # TODO: Maybe switch back to letting Ginga handle conversion
            #       to pixel coordinates.
            # Convert to pixels here (instead of in Ginga) because conversion
            # in Ginga was reportedly very slow.
            coord_x, coord_y = image.wcs.wcs.all_world2pix(
                coord_val.ra.deg, coord_val.dec.deg, 0)
            # In the event a *single* marker has been added, coord_x and coord_y
            # will be scalars. Make them arrays always.
            if np.ndim(coord_x) == 0:
                coord_x = np.array([coord_x])
                coord_y = np.array([coord_y])
        else:  # Use X,Y
            coord_x = table[x_colname].data
            coord_y = table[y_colname].data

        # Prepare canvas and retain existing marks
        try:
            c_mark = self._viewer.canvas.get_object_by_tag(marker_name)
        except Exception:
            objs = []
        else:
            objs = c_mark.objects
            self._viewer.canvas.delete_object_by_tag(marker_name)

        # TODO: Test to see if we can mix WCS and data on the same canvas
        objs += [self._marker(x=x, y=y, coord=coord_type)
                 for x, y in zip(coord_x, coord_y)]
        self._viewer.canvas.add(self.dc.CompoundObject(*objs), tag=marker_name)

    def remove_markers_by_name(self, marker_name):
        # TODO:
        #   arr : ``SkyCoord`` or array-like
        #   Sky coordinates or 2xN array.
        #
        # NOTE: How to match? Use np.isclose?
        #       What if there are 1-to-many matches?

        self.validate_marker_name(marker_name)
        if marker_name not in self._marktags:
            raise ValueError(
                f'Marker name {marker_name} not found in current markers. '
                f'Markers currently in use are {self.get_marker_names()}.')

        try:
            self._viewer.canvas.delete_object_by_tag(marker_name)
        except KeyError:
            self.logger.error(f'Unable to remove markers named {marker_name} from image.')
        else:
            self._marktags.remove(marker_name)

    def remove_all_markers(self):
        """Delete all markers using :meth:`remove_markers_by_name`."""
        # Grab the entire list of marker names before iterating
        # otherwise what we are iterating over changes.
        for marker_name in self.get_marker_names():
            self.remove_markers_by_name(marker_name)

    def validate_marker_name(self, marker_name):
        """Validate a given marker name.

        Parameters
        ----------
        marker_name : str
            Marker name to validate.

        Raises
        ------
        ValueError
            It is not allowed because the name is already defined in the
            ``RESERVED_MARKER_SET_NAMES`` attribute.

        """
        if marker_name in self.RESERVED_MARKER_SET_NAMES:
            raise ValueError(
                f"The marker name {marker_name} is not allowed. Any name is "
                f"allowed except these: {', '.join(self.RESERVED_MARKER_SET_NAMES)}")

    def set_cached_state(self):
        """Cache the following attributes before modifying their states:

        * ``click_center``
        * ``click_drag``
        * ``scroll_pan``

        This is used in :meth:`start_marking`, for example.
        """
        self._cached_state = dict(click_center=self.click_center,
                                  click_drag=self.click_drag,
                                  scroll_pan=self.scroll_pan)

    def restore_and_clear_cached_state(self):
        """Restore the following attributes with their cached states:

        * ``click_center``
        * ``click_drag``
        * ``scroll_pan``

        Then, clear the cache. This is used in :meth:`stop_marking`, for example.
        """
        self.click_center = self._cached_state['click_center']
        self.click_drag = self._cached_state['click_drag']
        self.scroll_pan = self._cached_state['scroll_pan']
        self._cached_state = {}

    @property
    def stretch_options(self):
        return self._viewer.get_color_algorithms()

    @property
    def stretch(self):
        """
        The image stretching algorithm in use.
        """
        return self._viewer.rgbmap.get_dist()

    # TODO: Possible to use astropy.visualization directly?
    @stretch.setter
    def stretch(self, value):
        valid_vals = self.stretch_options
        if value not in valid_vals:
            raise ValueError(f'Value must be one of: {valid_vals}')
        self._viewer.set_color_algorithm(value)

    @property
    def autocut_options(self):
        return self._viewer.get_autocut_methods()

    @property
    def cuts(self):
        return self._viewer.get_cut_levels()

    # TODO: Possible to use astropy.visualization directly?
    @cuts.setter
    def cuts(self, value):
        if isinstance(value, str):  # Autocut
            valid_vals = self.autocut_options
            if value not in valid_vals:
                raise ValueError(f'Value must be one of: {valid_vals}')
            self._viewer.set_autocut_params(value)
        else:
            if len(value) != 2:
                raise ValueError('Cut levels must be given as (low, high)')
            self._viewer.cut_levels(*value)

    @property
    def colormap_options(self):
        from ginga import cmap
        return cmap.get_names()

    def set_colormap(self, cmap):
        self._viewer.set_color_map(cmap)

    @property
    def cursor(self):
        """Current cursor information panel placement.

        Information must include the following:

        * X and Y cursor positions, depending on `pixel_offset`.
        * RA and Dec sky coordinates in HMS-DMS format, if available.
        * Value of the image under the cursor.

        You can set it to one of the following:

        * ``'top'`` places it above the image display.
        * ``'bottom'`` places it below the image display.
        * `None` hides it.

        """
        return self._cursor

    # NOTE: Subclass must re-implement if self._jup_coord is not ipyw.HTML
    #       or if self.ALLOWED_CURSOR_LOCATIONS is customized.
    @cursor.setter
    def cursor(self, value):
        if value is None:
            self._jup_coord.layout.visibility = 'hidden'
            self._jup_coord.layout.display = 'none'
        elif value in ('top', 'bottom'):
            self._jup_coord.layout.visibility = 'visible'
            self._jup_coord.layout.display = 'flex'
            if value == 'top':
                self.layout.flex_flow = 'column-reverse'
            else:
                self.layout.flex_flow = 'column'
        else:
            raise ValueError(
                f'Invalid value {value} for cursor. '
                f'Valid values are: {self.ALLOWED_CURSOR_LOCATIONS}')
        self._cursor = value

    @property

    def click_center(self):
        """When `True`, mouse left-click can be used to center an image.
        Otherwise, that interaction is disabled.

        You can set this property to `True` or `False`.
        This cannot be set to `True` when `is_marking` is also `True`.
        Setting this to `True` also disables `click_drag`.

        .. note:: In the future, this might accept non-bool values but not currently.

        """
        return self._click_center

    @click_center.setter
    def click_center(self, value):
        if not isinstance(value, bool):
            raise ValueError('Must be True or False')
        elif self.is_marking and value:
            raise ValueError('Interactive marking is in progress. Call '
                             'stop_marking() to end marking before setting '
                             'click_center')
        if value:
            self.click_drag = False

        self._click_center = value

    # Need this here because we need to overwrite the setter.
    @property
    def click_drag(self):
        """When `True`, the "click-and-drag" mode is an available interaction
        for panning. Otherwise, that interaction is disabled.

        You can set this property to `True` or `False`.
        This cannot be set to `True` when `is_marking` is also `True`.
        Setting this to `True` also disables `click_center`.

        """
        return self._click_drag

    @click_drag.setter
    def click_drag(self, value):
        if not isinstance(value, bool):
            raise ValueError('click_drag must be either True or False')
        if self.is_marking and value:
            raise ValueError('Interactive marking is in progress. Call '
                             'stop_marking() to end marking before setting '
                             'click_drag')
        self._click_drag = value
        bindmap = self._viewer.get_bindmap()
        if value:
            # Only turn off click_center if click_drag is being set to True
            self.click_center = False
            bindmap.map_event(None, (), 'ms_left', 'pan')
        else:
            bindmap.map_event(None, (), 'ms_left', 'cursor')

    # Need this here because we need to overwrite the setter.
    @property
    def scroll_pan(self):
        """When `True`, scrolling moves around (pans up/down) in the image.
        Otherwise, that interaction is disabled and becomes zoom.

        You can set this property to `True` or `False`.

        """
        return self._scroll_pan

    @scroll_pan.setter
    def scroll_pan(self, value):
        if not isinstance(value, bool):
            raise ValueError('scroll_pan must be either True or False')

        bindmap = self._viewer.get_bindmap()
        self._scroll_pan = value
        if value:
            bindmap.map_event(None, (), 'pa_pan', 'pan')
        else:
            bindmap.map_event(None, (), 'pa_pan', 'zoom')

    def save(self, filename, overwrite=False):
        if os.path.exists(filename) and not overwrite:
            raise ValueError(f'{filename} exists, use overwrite=True to force overwrite')

        ext = os.path.splitext(filename)[1].lower()
        if ext != '.png':
            raise ValueError(f'Extension {ext} not supported, use .png')

        # It turns out the image value is already in PNG format so we just
        # to write that out to a file.
        with open(filename, 'wb') as f:
            f.write(self._jup_img.value)


def _offset_is_pixel_or_sky(x):
    if isinstance(x, u.Quantity):
        if x.unit in (u.dimensionless_unscaled, u.pix):
            coord = 'data'
            val = x.value
        else:
            coord = 'wcs'
            val = x.to_value(u.deg)
    else:
        coord = 'data'
        val = x

    return val, coord
