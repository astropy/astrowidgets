from typing import Protocol, runtime_checkable, Any
from abc import abstractmethod
import os

from astropy.coordinates import SkyCoord
from astropy.nddata import NDData
from astropy.table import Table
from astropy.units import Quantity

from numpy.typing import ArrayLike

# Allowed locations for cursor display
ALLOWED_CURSOR_LOCATIONS = ('top', 'bottom', None)

# List of marker names that are for internal use only
RESERVED_MARKER_SET_NAMES = ('all',)

__all__ = [
    'ImageViewerInterface',
    'ALLOWED_CURSOR_LOCATIONS',
    'RESERVED_MARKER_SET_NAMES'
]


@runtime_checkable
class ImageViewerInterface(Protocol):
    # These are attributes, not methods. The type annotations are there
    # to make sure Protocol knows they are attributes. Python does not
    # do any checking at all of these types.
    click_center: bool
    click_drag: bool
    scroll_pan: bool
    image_width: int
    image_height: int
    zoom_level: float
    is_marking: bool
    stretch_options: tuple
    autocut_options: tuple
    cursor: str
    marker: Any
    cuts: Any
    stretch: str
    # viewer: Any

    # Allowed locations for cursor display
    ALLOWED_CURSOR_LOCATIONS: tuple = ALLOWED_CURSOR_LOCATIONS

    # List of marker names that are for internal use only
    RESERVED_MARKER_SET_NAMES: tuple = RESERVED_MARKER_SET_NAMES

    # The methods, grouped loosely by purpose

    # Methods for loading data
    @abstractmethod
    def load_fits(self, file: str | os.PathLike) -> None:
        """
        Load a FITS file into the viewer.

        Parameters
        ----------
        file : str or `astropy.io.fits.HDU`
            The FITS file to load. If a string, it can be a URL or a
            file path.
        """
        raise NotImplementedError

    @abstractmethod
    def load_array(self, array: ArrayLike) -> None:
        """
        Load a 2D array into the viewer.

        Parameters
        ----------
        array : array-like
            The array to load.
        """
        raise NotImplementedError

    @abstractmethod
    def load_nddata(self, data: NDData) -> None:
        """
        Load an `astropy.nddata.NDData` object into the viewer.

        Parameters
        ----------
        data : `astropy.nddata.NDData`
            The NDData object to load.
        """
        raise NotImplementedError

    # Saving contents of the view and accessing the view
    @abstractmethod
    def save(self, filename: str | os.PathLike, overwrite: bool = False) -> None:
        """
        Save the current view to a file.

        Parameters
        ----------
        filename : str
            The file to save to. The format is determined by the
            extension.

        overwrite : bool, optional
            If `True`, overwrite the file if it exists. Default is
            `False`.
        """
        raise NotImplementedError

    # Marker-related methods
    @abstractmethod
    def start_marking(self, marker_name: str | None = None) -> None:
        """
        Start interactive marking of points on the image.

        Parameters
        ----------
        marker_name : str, optional
            The name of the marker set to use. If not given, a unique
            name will be generated.
        """
        raise NotImplementedError

    @abstractmethod
    def stop_marking(self, clear_markers: bool = False) -> None:
        """
        Stop interactive marking of points on the image.

        Parameters
        ----------
        clear_markers : bool, optional
            If `True`, clear the markers that were created during
            interactive marking. Default is `False`.
        """
        raise NotImplementedError

    @abstractmethod
    def add_markers(self, table: Table, x_colname: str = 'x', y_colname: str = 'y',
                    skycoord_colname: str = 'coord', use_skycoord: bool = False,
                    marker_name: str | None = None) -> None:
        """
        Add markers to the image.

        Parameters
        ----------
        table : `astropy.table.Table`
            The table containing the marker positions.
        x_colname : str, optional
            The name of the column containing the x positions. Default
            is ``'x'``.
        y_colname : str, optional
            The name of the column containing the y positions. Default
            is ``'y'``.
        skycoord_colname : str, optional
            The name of the column containing the sky coordinates. If
            given, the ``use_skycoord`` parameter is ignored. Default
            is ``'coord'``.
        use_skycoord : bool, optional
            If `True`, the ``skycoord_colname`` column will be used to
            get the marker positions. Default is `False`.
        marker_name : str, optional
            The name of the marker set to use. If not given, a unique
            name will be generated.
        """
        raise NotImplementedError

    # @abstractmethod
    # def remove_all_markers(self):
    #     raise NotImplementedError

    @abstractmethod
    def reset_markers(self) -> None:
        """
        Remove all markers from the image.
        """
        raise NotImplementedError

    # @abstractmethod
    # def remove_markers_by_name(self, marker_name=None):
    #     raise NotImplementedError

    @abstractmethod
    def remove_markers(self, marker_name: str | None = None) -> None:
        """
        Remove markers from the image.

        Parameters
        ----------
        marker_name : str, optional
            The name of the marker set to remove. If not given, all
            markers will be removed.
        """
        raise NotImplementedError

    # @abstractmethod
    # def get_all_markers(self):
    #     raise NotImplementedError

    @abstractmethod
    def get_markers(self, x_colname: str = 'x', y_colname: str = 'y',
                    skycoord_colname: str = 'coord',
                    marker_name: str | None = None) -> Table:
        """
        Get the marker positions.

        Parameters
        ----------
        x_colname : str, optional
            The name of the column containing the x positions. Default
            is ``'x'``.
        y_colname : str, optional
            The name of the column containing the y positions. Default
            is ``'y'``.
        skycoord_colname : str, optional
            The name of the column containing the sky coordinates. Default
            is ``'coord'``.
        marker_name : str, optional
            The name of the marker set to use. If not given, all
            markers will be returned.

        Returns
        -------
        table : `astropy.table.Table`
            The table containing the marker positions.
        """
        raise NotImplementedError

    # Methods that modify the view
    @abstractmethod
    def center_on(self, point: tuple | SkyCoord):
        """
        Center the view on the point.

        Parameters
        ----------
        tuple or `~astropy.coordinates.SkyCoord`
            If tuple of ``(X, Y)`` is given, it is assumed
            to be in data coordinates.
        """
        raise NotImplementedError

    @abstractmethod
    def offset_by(self, dx: float | Quantity, dy: float | Quantity) -> None:
        """
        Move the center to a point that is given offset
        away from the current center.

        Parameters
        ----------
        dx, dy : float or `~astropy.units.Quantity`
            Offset value. Without a unit, assumed to be pixel offsets.
            If a unit is attached, offset by pixel or sky is assumed from
            the unit.
        """
        raise NotImplementedError

    @abstractmethod
    def zoom(self) -> None:
        """
        Zoom in or out by the given factor.

        Parameters
        ----------
        val : int
            The zoom level to zoom the image.
            See `zoom_level`.
        """
        raise NotImplementedError
