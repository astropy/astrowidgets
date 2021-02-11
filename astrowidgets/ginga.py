"""Module containing ``astrowidgets`` implemented with Ginga backend.

For this to work, ``astrowidgets`` must be installed along with the optional
dependencies specified for the Ginga backend; e.g.,::

    pip install 'astrowidgets[ginga]'

"""
import functools
import os
import warnings

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table

from ginga.AstroImage import AstroImage
from ginga.canvas.CanvasObject import drawCatalog
from ginga.web.jupyterw.ImageViewJpw import EnhancedCanvasView
from ginga.util.wcs import raDegToString, decDegToString

from astrowidgets.core import BaseImageWidget

__all__ = ['ImageWidget']


class ImageWidget(BaseImageWidget):
    """Image widget for Jupyter notebook using Ginga viewer.

    .. todo:: Any property passed to constructor has to be valid keyword.

    Parameters
    ----------
    logger : obj
        Ginga logger. For example::

            from ginga.misc.log import get_logger
            logger = get_logger('my_viewer', log_stderr=False,
                                log_file='ginga.log', level=40)

    use_opencv : bool, optional
        Let Ginga use ``opencv`` to speed up image transformation;
        e.g., rotation and mosaic. If this is enabled and you
        do not have ``opencv``, you will get a warning.

    image_height : int, optional
        Height of Jupyter notebook's image widget.

    image_width, pixel_coords_offset : int, optional
        See `~astrowidgets.core.BaseImageWidget`.

    """

    def __init__(self, logger, use_opencv=True, image_width=500,
                 image_height=500, pixel_coords_offset=0):
        super().__init__(image_width=image_width,
                         pixel_coords_offset=pixel_coords_offset)

        # TODO: Is this the best place for this?
        if use_opencv:
            try:
                from ginga import trcalc
                trcalc.use('opencv')
            except ImportError:
                logger.warning('install opencv or set use_opencv=False')

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
        return super().image_width

    @image_width.setter
    def image_width(self, value):
        # widgets expect width/height as strings, but most users will not, so
        # do the conversion.
        self._jup_img.width = str(value)
        self._viewer.set_window_size(self.image_width, self.image_height)

    # Need this here because we need to overwrite the setter.
    @property
    def image_height(self):
        return super().image_height

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

            val = (f'X: {data_x + self._pixel_offset:.2f}, '
                   f'Y: {data_y + self._pixel_offset:.2f}')

            if image.wcs.wcs is not None:
                try:
                    ra, dec = image.pixtoradec(data_x, data_y)
                    val += (f' (RA: {raDegToString(ra)},'
                            f' DEC: {decDegToString(dec)})')
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
            obj = self._marker(x=data_x, y=data_y)
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
                print(f'Centered on X={data_x + self._pixel_offset} '
                      f'Y={data_y + self._pixel_offset}')

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
            self._viewer.set_pan(*(np.asarray(point) - self._pixel_offset))

    def offset_to(self, dx, dy, skycoord_offset=False):
        if skycoord_offset:
            coord = 'wcs'
        else:
            coord = 'data'

        pan_x, pan_y = self._viewer.get_pan(coord=coord)
        self._viewer.set_pan(pan_x + dx, pan_y + dy, coord=coord)

    @property
    def zoom_level(self):
        return self._viewer.get_scale()

    zoom_level.__doc__ = (BaseImageWidget.zoom_level.__doc__
                          + "``'fit'`` means zoom to fit")

    @zoom_level.setter
    def zoom_level(self, value):
        if value == 'fit':
            self._viewer.zoom_fit()
        else:
            self._viewer.scale_to(value, value)

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
        return super().marker

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

        # Convert X,Y from 0-indexed to 1-indexed
        if self._pixel_offset != 0:
            xy_col += self._pixel_offset

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
            # Convert data coordinates from 1-indexed to 0-indexed
            if self._pixel_offset != 0:
                # Don't use the in-place operator -= here that modifies
                # the input table.
                coord_x = coord_x - self._pixel_offset
                coord_y = coord_y - self._pixel_offset

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

    @property
    def stretch_options(self):
        return self._viewer.get_color_algorithms()

    @property
    def stretch(self):
        return self._viewer.rgbmap.dist

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

    # Need this here because we need to overwrite the setter.
    @property
    def click_drag(self):
        return super().click_drag

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
        return super().scroll_pan

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

        self.logger.info(f'{filename} written')
