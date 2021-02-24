"""The ``astrowidgets.glupyter`` module contains a widget implemented with
the ``glue-jupyter`` (that uses ``bqplot``) backend.

For this to work, ``astrowidgets`` must be installed along with the optional
dependencies specified for the ``glue-jupyter`` (a.k.a. Glupyter or
``glupyter``) backend; e.g.,::

    pip install 'astrowidgets[glupyter]'

"""
import os
import warnings

import matplotlib.pyplot as plt
from astropy.io import fits

with warnings.catch_warnings():
    # Glue has warnings that we cannot do anything about.
    warnings.simplefilter('ignore')
    from bqplot_image_gl.interacts import MouseInteraction
    from glue.core.coordinates import (
        coordinates_from_header, coordinates_from_wcs)
    from glue.core.data import Data
    from glue_jupyter import jglue
    from glue_jupyter.bqplot.image.viewer import BqplotImageView
    from glue_jupyter.utils import validate_data_argument

from astrowidgets.core import BaseImageWidget

__all__ = ['ImageWidget']


class ImageWidget(BaseImageWidget):
    """Image widget for Jupyter notebook using ``glue-jupyter``/``bqplot``
    viewer.

    Parameters
    ----------
    show_axes : bool, optional
        If `True`, show X and Y axes.

    image_width : int, optional
        Width of Jupyter notebook's image widget.
        Height is automatically determined.

    pixel_coords_offset : int, optional
        An offset, typically either 0 or 1, to add/subtract to all
        pixel values when going to/from the displayed image.
        *In almost all situations the default value, ``0``, is the
        correct value to use.*

    """
    def __init__(self, show_axes=False, image_width=500, pixel_coords_offset=0):
        with warnings.catch_warnings():
            # Glue has warnings that we cannot do anything about.
            warnings.simplefilter('ignore')
            self._app = jglue()

        # Copied from glue.core.application_base
        self._viewer = BqplotImageView(self._app._session)
        self._viewer.register_to_hub(self._app._session.hub)
        self._app.add_widget(self._viewer)
        self._viewer.state.show_axes = show_axes

        # Cannot pass in viewer because ipywidgets want a basic widget as child
        super().__init__(image_widget=self._viewer.figure,
                         image_width=image_width,
                         pixel_coords_offset=pixel_coords_offset)

        # UNTIL HERE
        # TODO: This builds on PR 126. Example notebook in test_data/ztmp...

        # self.marker = ???

    @property
    def viewer(self):
        return self._jup_img

    @property
    def image_width(self):
        val = self.viewer.layout.width
        if val is None:
            return val
        try:
            v = int(val.strip('px'))
        except Exception:
            v = val  # Do not know how to convert to int
        return v

    @image_width.setter
    def image_width(self, value):
        if isinstance(value, (int, float)):
            value = f'{value}px'
        self.viewer.layout.width = value

    @property
    def image_height(self):
        val = self.viewer.layout.height
        if val is None:
            return val
        try:
            v = int(val.strip('px'))
        except Exception:
            v = val  # Do not know how to convert to int
        return v

    @image_height.setter
    def image_height(self, value):
        if isinstance(value, (int, float)):
            value = f'{value}px'
        self.viewer.layout.height = value

    def _link_image_to_cb(self):
        image = self.viewer.marks[0]
        interaction = MouseInteraction(
            x_scale=image.scales['x'], y_scale=image.scales['y'],
            move_throttle=70)
        self.viewer.interaction = interaction
        interaction.on_msg(self._mouse_cb)

    def _mouse_cb(self, interaction, data, buffers):
        """Callback for mouse events; e.g., to display position in RA/DEC deg."""
        event = data['event']

        if event not in ('mousemove', 'click'):
            return  # no-op

        # Needs https://github.com/glue-viz/bqplot-image-gl/pull/37 .
        # Logic based on Tom Robitaille's Jdaviz/Imviz demo.
        if event == 'mousemove':
            image = self._viewer.state.reference_data  # Get active image

            # Extract data coordinates - these are pixels in the image
            x = data['domain']['x']
            y = data['domain']['y']
            val = f'X: {x + self._pixel_offset:.2f}, Y: {y + self._pixel_offset:.2f}'

            # Convert these to a SkyCoord via WCS
            if image.coords is not None:
                sky = image.coords.pixel_to_world(x, y).icrs
                ra = sky.ra.to_string(sep='hms')
                dec = sky.dec.to_string(sep='hms')
                val += f' (RA: {ra}, DEC: {dec})'

            # Extract data values at this position
            ix = int(round(x))
            iy = int(round(y))
            if ix >= 0 and iy >= 0 and ix < image.shape[1] and iy < image.shape[0]:
                value = image.get_data(image.main_components[0])[iy, ix]
            else:
                value = "NA"
            val += f', value: {value}'

        # FIXME: Handle clicks
        elif event == 'click':
            x = data['domain']['x']
            y = data['domain']['y']

            if self.is_marking:
                marker_name = self._interactive_marker_set_name

                # FIXME: Marking logic here, how to display marks and how to retrieve them later?

                # For debugging.
                with self.print_out:
                    print(f'Selected {x} {y}')

            elif self.click_center:
                # FIXME: Implement this
                #self.center_on((x, y))

                # For debugging.
                with self.print_out:
                    print(f'Centered on X={x + self._pixel_offset} '
                          f'Y={y + self._pixel_offset}')

        self._jup_coord.value = val

    def load_fits(self, filename, numhdu=0, memmap=None):
        """Load a FITS file into the viewer.

        Parameters
        ----------
        filename : str
            Name of the FITS file.

        numhdu : int
            Extension number of the desired HDU. If not given, it is
            assumed to be in EXT 0.

        memmap : bool or `None`
            Memory mapping. See `astropy.io.fits.open`.

        """
        with fits.open(filename, memmap=memmap) as pf:
            data = _new_data_from_hdulist(pf, numhdu)

        viewer = self._viewer
        data = validate_data_argument(self._app.data_collection, data)

        # NOTE: Changing state via API does not work properly if we loaded
        # multiple data into the same widget in the same session.
        self._app.data_collection.clear()

        self._app.add_datasets(data)
        viewer.add_data(data)
        viewer.state.reference_data = data  # Bring it to focus
        self._link_image_to_cb()

    def load_nddata(self, nddata, label=None):
        """Load a `~astropy.nddata.NDData` object into the viewer.

        .. todo:: Add flag/masking support, etc.

        Parameters
        ----------
        nddata : `~astropy.nddata.NDData`
            ``NDData`` with image data and WCS.

        label : str or `None`
            Glue label name for the data.

        """
        if label is None:
            label = _gen_random_label()

        data = Data(label=label)

        if nddata.wcs is not None:
            coords = coordinates_from_wcs(nddata.wcs)
            data.coords = coords

        data.add_component(component=nddata.data, label=label)

        viewer = self._viewer
        data = validate_data_argument(self._app.data_collection, data)

        # NOTE: Changing state via API does not work properly if we loaded
        # multiple data into the same widget in the same session.
        self._app.data_collection.clear()

        self._app.add_datasets(data)
        viewer.add_data(data)
        viewer.state.reference_data = data  # Bring it to focus
        self._link_image_to_cb()

    def load_array(self, arr, label=None):
        """Load a 2D array into the viewer.

        .. note:: Use :meth:`load_nddata` for WCS support.

        Parameters
        ----------
        arr : array-like
            2D array.

        label : str or `None`
            Glue label name for the data.

        """
        if label is None:
            label = _gen_random_label()

        data = Data(label=label)
        data.add_component(component=arr, label=label)

        viewer = self._viewer
        data = validate_data_argument(self._app.data_collection, data)

        # NOTE: Changing state via API does not work properly if we loaded
        # multiple data into the same widget in the same session.
        self._app.data_collection.clear()

        self._app.add_datasets(data)
        viewer.add_data(data)
        viewer.state.reference_data = data  # Bring it to focus
        self._link_image_to_cb()

    def center_on(self, point):
        raise NotImplementedError  # FIXME

    def offset_to(self, dx, dy, skycoord_offset=False):
        raise NotImplementedError  # FIXME

    @property
    def zoom_level(self):
        raise NotImplementedError  # FIXME

    @zoom_level.setter
    def zoom_level(self, value):
        raise NotImplementedError  # FIXME

    # Need this here because we need to overwrite the setter.
    @property
    def marker(self):
        return super().marker

    @marker.setter
    def marker(self, value):
        raise NotImplementedError  # FIXME

    def get_markers_by_name(self, marker_name, x_colname='x', y_colname='y',
                            skycoord_colname='coord'):
        raise NotImplementedError  # FIXME

    def add_markers(self, table, x_colname='x', y_colname='y',
                    skycoord_colname='coord', use_skycoord=False,
                    marker_name=None):
        raise NotImplementedError  # FIXME

    def remove_markers_by_name(self, marker_name):
        raise NotImplementedError  # FIXME

    @property
    def stretch_options(self):
        # Options extracted from Glue
        return ['arcsinh', 'linear', 'log', 'sqrt']

    @property
    def stretch(self):
        return self._viewer.state.layers[0].stretch

    @stretch.setter
    def stretch(self, value):
        self._viewer.state.layers[0].stretch = value

    @property
    def autocut_options(self):
        raise NotImplementedError  # FIXME: Does Glue support this?

    @property
    def cuts(self):
        return (self._viewer.state.layers[0].v_min,
                self._viewer.state.layers[0].v_max)

    @cuts.setter
    def cuts(self, value):
        if len(value) != 2:
            raise ValueError('Cut levels must be given as (low, high)')
        self._viewer.state.layers[0].v_min = value[0]
        self._viewer.state.layers[0].v_max = value[1]

    @property
    def colormap_options(self):
        return plt.colormaps()

    def set_colormap(self, cmap):
        self._viewer.state.layers[0].cmap = plt.cm.get_cmap(cmap)

    def save(self, filename, **kwargs):
        """Save the current image view to given filename (PNG only).

        Parameters
        ----------
        filename : str
            Image filename. A window dialog will pop up for directory
            selection.

        kwargs : dict, optional
            Not used.

        Raises
        ------
        ValueError
            Invalid input.

        """
        ext = os.path.splitext(filename)[1].lower()
        if ext != '.png':
            raise ValueError(f'Extension {ext} not supported, use .png')

        self.viewer.save_png(filename)


def _new_data_from_hdulist(hdulist, numhdu):
    """This code is extracted out of Glue fits_reader."""
    hdulist_name = hdulist.filename()
    if hdulist_name is None:
        hdulist_name = 'HDUList'

    label_base = os.path.basename(hdulist_name).rpartition('.')[0]

    if not label_base:
        label_base = os.path.basename(hdulist_name)

    label = f'{label_base}[{numhdu}]'
    data = Data(label=label)
    hdu = hdulist[numhdu]
    coords = coordinates_from_header(hdu.header)
    data.coords = coords
    data.add_component(component=hdu.data, label=numhdu)
    return data


def _gen_random_label(prefix='data|'):
    """Generate unique identifier."""
    import base64
    import uuid

    return f'{prefix}{str(base64.b85encode(uuid.uuid4().bytes), "utf-8")}'
