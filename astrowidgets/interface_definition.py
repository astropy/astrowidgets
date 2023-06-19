from typing import Protocol, runtime_checkable, Any
from abc import abstractmethod

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
    # This are attributes, not methods. The type annotations are there
    # to make sure Protocol knows they are attributes. Python does not
    # do any checking at all of these types.
    click_center: bool
    click_drag: bool
    scroll_pan: bool
    image_width: int
    image_height: int
    zoom_level: float
    marker: Any
    cuts: Any
    stretch: str
    # viewer: Any

    # Allowed locations for cursor display
    ALLOWED_CURSOR_LOCATIONS : tuple = ALLOWED_CURSOR_LOCATIONS

    # List of marker names that are for internal use only
    RESERVED_MARKER_SET_NAMES : tuple = RESERVED_MARKER_SET_NAMES

    # The methods, grouped loosely by purpose

    # Methods for loading data
    @abstractmethod
    def load_fits(self, file):
        raise NotImplementedError

    @abstractmethod
    def load_array(self, array):
        raise NotImplementedError

    @abstractmethod
    def load_nddata(self, data):
        raise NotImplementedError

    # Saving contents of the view and accessing the view
    @abstractmethod
    def save(self, filename):
        raise NotImplementedError

    # Marker-related methods
    @abstractmethod
    def start_marking(self, marker_name=None):
        raise NotImplementedError

    @abstractmethod
    def stop_marking(self, clear_markers=False):
        raise NotImplementedError

    @abstractmethod
    def add_markers(self, table, x_colname='x', y_colname='y',
                    skycoord_colname='coord', use_skycoord=False,
                    marker_name=None):
        raise NotImplementedError

    # @abstractmethod
    # def remove_all_markers(self):
    #     raise NotImplementedError

    @abstractmethod
    def reset_markers(self):
        raise NotImplementedError

    # @abstractmethod
    # def remove_markers_by_name(self, marker_name=None):
    #     raise NotImplementedError

    @abstractmethod
    def remove_markers(self, marker_name=None):
        raise NotImplementedError

    # @abstractmethod
    # def get_all_markers(self):
    #     raise NotImplementedError

    @abstractmethod
    def get_markers(self, x_colname='x', y_colname='y',
                    skycoord_colname='coord',
                    marker_name=None):
        raise NotImplementedError

    # Methods that modify the view
    @abstractmethod
    def center_on(self):
        raise NotImplementedError

    @abstractmethod
    def offset_by(self):
        raise NotImplementedError

    @abstractmethod
    def zoom(self):
        raise NotImplementedError
