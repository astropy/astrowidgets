from typing import Protocol, runtime_checkable, Any
from abc import abstractmethod


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
    def start_marking(self):
        raise NotImplementedError

    @abstractmethod
    def stop_marking(self):
        raise NotImplementedError

    @abstractmethod
    def add_markers(self):
        raise NotImplementedError

    @abstractmethod
    def remove_all_markers(self):
        raise NotImplementedError

    @abstractmethod
    def remove_markers_by_name(self, marker_name=None):
        raise NotImplementedError

    @abstractmethod
    def get_all_markers(self):
        raise NotImplementedError

    @abstractmethod
    def get_markers_by_name(self, marker_name=None):
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
