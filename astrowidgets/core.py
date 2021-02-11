"""Module containing an abstract base class defining the core functionality
of ``astrowidgets``.

A working implementation using Ginga is in `astrowidgets.ginga.ImageWidget`.

"""
from abc import ABCMeta, abstractmethod

from astropy.coordinates import SkyCoord
from astropy.table import vstack

import ipywidgets as ipyw

__all__ = ['BaseImageWidget']


class _MetaWidget(ABCMeta, type(ipyw.VBox)):
    pass


class BaseImageWidget(ipyw.VBox, metaclass=_MetaWidget):
    """Base class for image widget for Jupyter notebook.

    .. note::

        The constructor of subclass must assign ``self.marker`` after
        calling ``super().__init__(...)``; also see `marker`.

        Subclass must implement necessary mouse and keyboard event
        handlings not laid out here to make some of the documented methods,
        attributes, and properties (e.g., `click_center`) work.

        In some cases, subclass might also want to re-implement non-abstract
        API as needed, using the implementation in this base class as a guide.

    Parameters
    ----------
    image_widget : obj or `None`, optional
        Widget for image display. If not given, ``ipywidgets.Image`` is used.

    cursor_widget : obj or `None`, optional
        Widget for cursor information display. If not given,
        ``ipywidgets.HTML`` is used.

    image_width : int, optional
        Width of Jupyter notebook's image widget.
        Height is automatically determined.

    pixel_coords_offset : int, optional
        An offset, typically either 0 or 1, to add/subtract to all
        pixel values when going to/from the displayed image.
        *In almost all situations the default value, ``0``, is the
        correct value to use.*

    Attributes
    ----------
    children : list
        Children of the ``ipywidgets.VBox`` that consists of the given
        ``image_widget`` and ``cursor_widget`` to be displayed.

    print_out : obj
        ``ipywidgets.Output`` instance for printing output to a notebook
        cell; this is especially useful for debugging.

    ALLOWED_CURSOR_LOCATIONS : list
        Possible `cursor` widget placements relative to image widget.

    RESERVED_MARKER_SET_NAMES : list
        Marker names reserved for internal use only.

    """
    def __init__(self, image_widget=None, cursor_widget=None, image_width=500,
                 pixel_coords_offset=0):
        super().__init__()

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

        self._pixel_offset = pixel_coords_offset

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

    @property
    @abstractmethod
    def viewer(self):
        """The underlying viewer tied to the backend."""
        pass

    @property
    def image_width(self):
        """Width of image widget."""
        return int(self._jup_img.width)

    @image_width.setter
    def image_width(self, value):
        # widgets expect width/height as strings, but most users will not, so
        # do the conversion.
        self._jup_img.width = str(value)

    @property
    def image_height(self):
        """Height of image widget."""
        return int(self._jup_img.height)

    @image_height.setter
    def image_height(self, value):
        # widgets expect width/height as strings, but most users will not, so
        # do the conversion.
        self._jup_img.height = str(value)

    @property
    def pixel_offset(self):
        """An offset, typically either 0 or 1, to add/subtract to all
        pixel values when going to/from the displayed image.
        *In almost all situations the default value, ``0``, is the
        correct value to use.*

        This value cannot be modified after initialization.

        """
        return self._pixel_offset

    @abstractmethod
    def load_fits(self, filename, **kwargs):
        """Load a FITS file into the viewer.

        Parameters
        ----------
        filename : str
            Name of the FITS file.

        kwargs : dict, optional
            Keywords for the loader specific to the chosen backend, if any.

        """
        pass

    @abstractmethod
    def load_nddata(self, nddata, **kwargs):
        """Load a `~astropy.nddata.NDData` object into the viewer.

        Parameters
        ----------
        nddata : `~astropy.nddata.NDData`
            ``NDData`` with image data and WCS.

        kwargs : dict, optional
            Keywords for the loader specific to the chosen backend, if any.

        """
        pass

    @abstractmethod
    def load_array(self, arr, **kwargs):
        """Load a 2D array into the viewer.

        .. note:: Use :meth:`load_nddata` for WCS support.

        Parameters
        ----------
        arr : array-like
            2D array.

        kwargs : dict, optional
            Keywords for the loader specific to the chosen backend, if any.

        """
        pass

    @abstractmethod
    def center_on(self, point):
        """Centers the view on a particular point.

        Parameters
        ----------
        point : tuple or `~astropy.coordinates.SkyCoord`
            If tuple of ``(X, Y)`` is given, it is assumed
            to be in data coordinates. If data coordinates is given,
            `pixel_offset` needs to be taken into account as well.

        """
        pass

    @abstractmethod
    def offset_to(self, dx, dy, skycoord_offset=False):
        """Move the center to a point that is given offset
        away from the current center.

        Parameters
        ----------
        dx, dy : float
            Offset values. Unit is assumed based on
            ``skycoord_offset``.

        skycoord_offset : bool, optional
            If `True`, offset must be given in degrees.
            Otherwise, they are in pixel values.

        """
        pass

    @property
    @abstractmethod
    def zoom_level(self):
        """Zoom level (settable):

        * 1 means real-pixel-size.
        * 2 means zoomed in by a factor of 2.
        * 0.5 means zoomed out by a factor of 2.

        """
        pass

    @zoom_level.setter
    @abstractmethod
    def zoom_level(self, value):
        pass

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

        .. note::

            Subclass should overwrite this docstring with examples specific to its backend; e.g., using ``marker.__doc__``.

        .. note::

            Its setter must set the following:

            * ``self._marker_dict`` with the given dictionary.
            * ``self._marker`` with an **object** built from the given dictionaru for the backend.

        """
        return self._marker_dict

    @marker.setter
    @abstractmethod
    def marker(self, value):
        # See notes for the marker property.
        pass

    def get_marker_names(self):
        """Return a list of used marker names.

        Returns
        -------
        names : list of str
            Sorted list of marker names.

        """
        return sorted(self._marktags)

    @abstractmethod
    def get_markers_by_name(self, marker_name, x_colname='x', y_colname='y',
                            skycoord_colname='coord'):
        """Return the locations of markers for the given name.

        Parameters
        ----------
        marker_name : str
            Available names can be obtained by calling :meth:`get_marker_names`.

        x_colname, y_colname : str, optional
            Column names for X and Y data coordinates.
            Coordinates returned are 0- or 1-indexed, depending
            on `pixel_offset`.

        skycoord_colname : str, optional
            Column name for `~astropy.coordinates.SkyCoord`, which contains
            sky coordinates associated with the active image.
            This is ignored if image has no WCS.

        Returns
        -------
        markers_table : `~astropy.table.Table` or `None`
            Table of markers, if any, contains the following columns:

            * x (or as set by ``x_colname``)
            * y (or as set by ``y_colname``)
            * (OPTIONAL) coord (or as set by ``skycoord_colname``) -- Only if available
            * marker name (from ``marker_name``) -- Useful for :meth:`get_all_markers`

        Raises
        ------
        ValueError
            Marker name is invalid.

        See also
        --------
        get_all_markers

        """
        pass

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

    @abstractmethod
    def add_markers(self, table, x_colname='x', y_colname='y',
                    skycoord_colname='coord', use_skycoord=False,
                    marker_name=None):
        """Show markers in the image at given points using the current
        marker style.

        Parameters
        ----------
        table : `~astropy.table.Table`
            Table containing marker locations. Compulsory columns depend on
            ``use_skycoord``.

        x_colname, y_colname : str, optional
            Column names for X and Y. Coordinates can be 0- or 1-indexed, as
            given by `pixel_offset`. These are only used if
            ``use_skycoord=False``.

        skycoord_colname : str, optional
            Column name with `~astropy.coordinates.SkyCoord` objects.
            This is only used if ``use_skycoord=True``.

        use_skycoord : bool, optional
            If `False`, mark using ``x_colname`` and ``y_colname``;
            otherwise ``skycoord_colname``.

        marker_name : str or `None`, optional
            Name to assign the markers in the table. If not given, an internal
            default is used (set by the ``_default_mark_tag_name`` attribute).
            A given name cannot fail :meth:`validate_marker_name` check.
            Marker name used will be added to ``_marktags`` attribute.

        Raises
        ------
        ValueError
            Marker name is invalid or operation failed.

        See also
        --------
        get_markers_by_name

        """
        pass

    @abstractmethod
    def remove_markers_by_name(self, marker_name):
        """Remove all of the markers by the name used on addition.

        Parameters
        ----------
        marker_name : str
            Name used when the markers were added. Available names can be
            obtained by calling :meth:`get_marker_names`.
            This name will be removed from ``_marktags`` attribute.

        Raises
        ------
        ValueError
            Marker name is invalid.

        See also
        --------
        remove_all_markers

        """
        pass

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
    @abstractmethod
    def stretch_options(self):
        """List of all available options for image stretching."""
        pass

    @property
    @abstractmethod
    def stretch(self):
        """The image stretching algorithm in use."""
        pass

    @stretch.setter
    @abstractmethod
    def stretch(self, value):
        pass

    @property
    @abstractmethod
    def autocut_options(self):
        """List of all available options for image auto-cut."""
        pass

    @property
    @abstractmethod
    def cuts(self):
        """Current image cut levels as ``(low, high)``.

        To set new cut levels, provide one of the following:

        * A tuple of ``(low, high)`` values.
        * One of the options returned by `autocut_options`.

        """
        pass

    @cuts.setter
    @abstractmethod
    def cuts(self, value):
        pass

    @property
    @abstractmethod
    def colormap_options(self):
        """List of available colormap names."""
        pass

    @abstractmethod
    def set_colormap(self, cmap):
        """Set colormap to the given name.

        Parameters
        ----------
        cmap : str
            Colormap name. Possible values can be obtained from
            :meth:`colormap_options`.

        """
        pass

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
        if value:
            # Only turn off click_center if click_drag is being set to True
            self.click_center = False

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
        self._scroll_pan = value

    @abstractmethod
    def save(self, filename, overwrite=False):
        """Save the current image view to given filename.
        File type (e.g., PNG) is assumed from the given extension.

        Parameters
        ----------
        filename : str
            Image filename. If you want to save it in a different
            directory, provide the full path.

        overwrite : bool, optional
            Overwrite existing file with the same name.

        Raises
        ------
        ValueError
            Invalid input.

        """
        pass
