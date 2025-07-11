import inspect
from pathlib import Path

import numpy as np
from astropy.nddata import NDData
import astropy.visualization as apviz

from bqplot import Figure, LinearScale, Axis, ColorScale, PanZoom, ScatterGL
from bqplot_image_gl import ImageGL
from bqplot_image_gl.interacts import (MouseInteraction,
                                       keyboard_events, mouse_events)

import ipywidgets as ipw

from matplotlib import pyplot
from matplotlib.colors import to_hex

from astro_image_display_api import ImageViewerInterface
from astro_image_display_api.image_viewer_logic import ImageViewerLogic

import traitlets as trait


def docs_from_super_if_missing(cls):
    """
    Decorator to copy the docstrings from the interface methods to the
    methods in the class.
    """
    for name, method in cls.__dict__.items():
        if not name.startswith("_"):
            if method.__doc__:
                continue
            # Not sure why this fails, but it does.
            # method.__doc__ = inspect.getdoc(method)
            interface_method = getattr(ImageViewerLogic, name, None)

            if interface_method:
                method.__doc__ = interface_method.__doc__
            # print(f"{method} {method.__doc__} {inspect.getdoc(method)} {interface_method.__doc__=}")
            # supers = cls.__mro__
            # for a_super in supers:
            #     print(f"{name=} {a_super=}")
            #     interface_method = getattr(a_super, name, None)
            #     if interface_method:
            #         print(f"{name=} {a_super=}")
            #         method.__doc__ = interface_method.__doc__
            #         break
    return cls


class _AstroImage(ipw.VBox):
    """
    Encapsulate an image as a bqplot figure inside a box.

    bqplot is involved for its pan/zoom capabilities, and it presents as
    a box to obscure the usual bqplot properties and methods.
    """
    def __init__(self, image_data=None,
                 display_width=500,
                 viewer_aspect_ratio=1.0):
        super().__init__()

        self._viewer_aspect_ratio = viewer_aspect_ratio

        self._display_width = display_width
        self._display_height = self._viewer_aspect_ratio * self._display_width


        layout = ipw.Layout(width=f'{self._display_width}px',
                            height=f'{self._display_height}px',
                            justify_content='center')

        self._figure_layout = layout

        scale_x = LinearScale(min=0, max=1, #self._image_shape[1],
                              allow_padding=False)
        scale_y = LinearScale(min=0, max=1, #self._image_shape[0],
                              allow_padding=False)
        self._scales = {'x': scale_x, 'y': scale_y}
        axis_x = Axis(scale=scale_x, visible=False)
        axis_y = Axis(scale=scale_y, orientation='vertical', visible=False)
        scales_image = {'x': scale_x, 'y': scale_y,
                        'image': ColorScale(max=1, min=0,
                                            scheme='Greys')}

        self._image_shape = None

        self._scatter_marks = {}

        self._figure = Figure(scales=self._scales, axes=[axis_x, axis_y],
                              fig_margin=dict(top=0, left=0,
                                              right=0, bottom=0),
                              layout=layout)

        self._image = ImageGL(image=np.zeros(shape=[1, 1]), scales=scales_image)

        self._figure.marks = (self._image, )

        panzoom = PanZoom(scales={'x': [scales_image['x']],
                                  'y': [scales_image['y']]})
        interaction = MouseInteraction(x_scale=scales_image['x'],
                                       y_scale=scales_image['y'],
                                       move_throttle=70, next=panzoom,
                                       events=keyboard_events + mouse_events)

        self._figure.interaction = interaction

        # Keep track of this separately so that it is easy to change
        # its state.
        self._panzoom = panzoom

        if image_data:
            self.set_data(image_data, reset_view=True)

        self.children = (self._figure, )

    @property
    def data_aspect_ratio(self):
        """
        Aspect ratio of the image data, horizontal size over vertical size.
        """
        return self._image_shape[0] / self._image_shape[1]

    def reset_scale_to_fit_image(self):
        wide = self.data_aspect_ratio < 1
        tall = self.data_aspect_ratio > 1
        square = self.data_aspect_ratio == 1

        if wide:
            self._scales['x'].min = 0
            self._scales['x'].max = self._image_shape[1]
            self._set_scale_aspect_ratio_to_match_viewer()
        elif tall or square:
            self._scales['y'].min = 0
            self._scales['y'].max = self._image_shape[0]
            self._set_scale_aspect_ratio_to_match_viewer(reset_scale='x')

        # Great, now let's center
        self.center = (self._image_shape[1]/2,
                       self._image_shape[0]/2)


    def _set_scale_aspect_ratio_to_match_viewer(self,
                                                reset_scale='y'):
        # Set the scales so that they match the aspect ratio
        # of the viewer, preserving the current image center.
        width_x, width_y = self.scale_widths
        frozen_width = dict(y=width_x, x=width_y)
        scale_aspect = width_x / width_y
        figure_x = float(self._figure.layout.width[:-2])
        figure_y = float(self._figure.layout.height[:-2])
        figure_aspect = figure_x / figure_y
        current_center = self.center
        if abs(figure_aspect - scale_aspect) > 1e-4:
            # Make the scale aspect ratio match the
            # figure layout aspect ratio
            if reset_scale == 'y':
                scale_factor = 1/ figure_aspect
            else:
                scale_factor = figure_aspect

            self._scales[reset_scale].min = 0
            self._scales[reset_scale].max = frozen_width[reset_scale] * scale_factor
            self.center = current_center

    def set_data(self, image_data, reset_view=True):
        self._image_shape = image_data.shape

        if reset_view:
            self.reset_scale_to_fit_image()

        # Set the image data and map it to the bqplot figure so that
        # cursor location corresponds to the underlying array index.
        # The offset follows the convention that the index corresponds
        # to the center of the pixel.
        self._image.image = image_data
        self._image.x = [-0.5, self._image_shape[1] - 0.5]
        self._image.y = [-0.5, self._image_shape[0] - 0.5]

    @property
    def scale_widths(self):
        width_x = self._scales['x'].max - self._scales['x'].min
        width_y = self._scales['y'].max - self._scales['y'].min
        return (width_x, width_y)

    @property
    def viewer_size(self):
        """
        The size of the viewer in pixels, as a tuple of (width, height).
        """
        width = float(self._figure.layout.width[:-2])
        height = float(self._figure.layout.height[:-2])
        return (width, height)

    @property
    def center(self):
        """
        Center of current view in pixels in x, y.
        """
        x_center = (self._scales['x'].min + self._scales['x'].max) / 2
        y_center = (self._scales['y'].min + self._scales['y'].max) / 2
        return (x_center, y_center)

    @center.setter
    def center(self, value):
        x_c, y_c = value

        width_x, width_y = self.scale_widths
        self._scales['x'].max = x_c + width_x / 2
        self._scales['x'].min = x_c - width_x / 2
        self._scales['y'].max = y_c + width_y / 2
        self._scales['y'].min = y_c - width_y / 2

    @property
    def interaction(self):
        return self._figure.interaction

    def set_color(self, colors):
        # colors here means a list of hex colors
        self._image.scales['image'].colors = colors

    def _check_file_exists(self, filename, overwrite=False):
        if Path(filename).exists() and not overwrite:
            raise ValueError(f'File named {filename} already exists. Use '
                             f'overwrite=True to overwrite it.')

    def save_png(self, filename, overwrite=False):
        self._check_file_exists(filename, overwrite=overwrite)
        self._figure.save_png(filename)

    def save_svg(self, filename, overwrite=False):
        self._check_file_exists(filename, overwrite=overwrite)
        self._figure.save_svg(filename)

    def set_pan(self, on_or_off):
        self._panzoom.allow_pan = on_or_off

    def set_scroll_zoom(self, on_or_off):
        self._panzoom.allow_zoom = on_or_off

    def set_size(self, size, direction):
        """
        Set the size of the scales to the desired size.

        Parameters
        ----------
        size : float
            The size in pixels to set the scale to.

        direction : {'x', 'y', 'smallest'}
            The direction in which to set the size. If 'x', the x scale
            will be set to the desired size. If 'y', the y scale will be
            set to the desired size. If 'smallest', the scale with the
            smallest viewer size will be set to the desired size.
        """
        if direction == "smallest":
            # Set the scale with the smallest width to the desired size
            if self.viewer_size[0] < self.viewer_size[1]:
                direction = 'x'
            else:
                direction = 'y'

        scale_to_set = self._scales[direction]
        cen = {}
        cen['x'], cen['y'] = self.center
        with scale_to_set.hold_trait_notifications():
            scale_to_set.min = cen[direction] - size/2
            scale_to_set.max = cen[direction] + size/2

        reset_scale = 'x' if direction == 'y' else 'y'

        self._set_scale_aspect_ratio_to_match_viewer(reset_scale)

    def set_zoom_level(self, zoom_level):
        """
        Set zoom level of viewer. A zoom level of 1 means 1 pixel
        in the image is 1 pixel in the viewer, i.e. the scale width
        in the horizontal direction matches the width in pixels
        of the figure.
        """

        # The width is reset here but the height could be set instead
        # and the result would be the same.
        figure_width = float(self._figure.layout.width[:-2])
        new_width = figure_width / zoom_level
        self.set_size(new_width, 'x')
        self._set_scale_aspect_ratio_to_match_viewer('y')

    def get_current_width(self):
        """
        Get the zoom level of the current view, if such a view has been set.

        A zoom level of 1 means 1 pixel in the image is 1 pixel in the viewer,
        i.e. the scale width in the horizontal direction matches the width in
        pixels of the figure.
        """
        if self._image_shape is None:
            return None

        # The width is used here but the height could be used instead
        # and the result would be the same since the pixels are square.
        scale_width = self.scale_widths[1]

        return scale_width

    def plot_named_markers(self, x, y, mark_id, color='yellow',
                           size=100, shape='circle', **kwd):
        scale_dict = dict(x=self._scales['x'], y=self._scales['y'])
        sc = ScatterGL(scales=scale_dict,
                       x=x, y=y,
                       colors=[color],
                       default_size=100,
                       marker=shape,
                       fill=False)

        self._scatter_marks[mark_id] = sc
        self._update_marks()

    def remove_named_markers(self, mark_id):
        if isinstance(mark_id, str):
            if mark_id == '*':
                self.remove_markers()
                return
            else:
                mark_id = [mark_id]
        for m_id in mark_id:
            try:
                del self._scatter_marks[m_id]
            except KeyError:
                raise ValueError(f'Markers {m_id} are not present.')

        self._update_marks()

    def remove_markers(self):
        self._scatter_marks = {}
        self._update_marks()

    def _update_marks(self):
        marks = [self._image] + [mark for mark in self._scatter_marks.values()]
        self._figure.marks = marks


def bqcolors(colormap, reverse=False):
    # bqplot-image-gl has 256 levels
    LEVELS = 256

    # Make a matplotlib colormap object
    mpl = pyplot.get_cmap(colormap, LEVELS)

    # Get RGBA colors
    mpl_colors = mpl(np.linspace(0, 1, LEVELS))

    # Convert RGBA to hex
    bq_colors = [to_hex(mpl_colors[i, :]) for i in range(LEVELS)]

    if reverse:
        bq_colors = bq_colors[::-1]

    return bq_colors


# The inheritance order below matters -- VBox needs to come first
@docs_from_super_if_missing
class ImageWidget(ipw.VBox, ImageViewerLogic):
    def __init__(self, *args, display_width=500, display_aspect_ratio=1):
        super().__init__(*args)
        self._set_up_catalog_image_dicts()
        # self.image_width = image_width
        # self.image_height = image_height

        self._astro_im = _AstroImage(display_width=display_width,
                                     viewer_aspect_ratio=display_aspect_ratio)
        self._default_cuts = apviz.MinMaxInterval()
        self._default_stretch = None

        self._data = None
        self._wcs = None

        # Use this to manage whether or not to send changes in zoom level
        # to the viewer.
        self._viewport_change_source_is_gui = False

        # Provide an Output widget to which prints can be directed for
        # debugging.
        self._print_out = ipw.Output()

        self.marker = {'color': 'red', 'radius': 20, 'type': 'square'}
        self._cuts = apviz.AsymmetricPercentileInterval(1, 99)

        self._cursor = ipw.HTML('Coordinates show up here')

        self._init_mouse_callbacks()
        self._init_watch_image_changes()
        self.children = [self._astro_im, self._cursor]

    def _init_mouse_callbacks(self):

        def on_mouse_message(interaction, event_data, buffers):
            """
            This function simply detects the event type then dispatches
            to the method that handles that event.

            The ``event_data`` contains all of the information we need.
            """
            if event_data['event'] == 'mousemove':
                self._mouse_move(event_data)
            elif event_data['event'] == 'click':
                self._mouse_click(event_data)

        self._astro_im.interaction.on_msg(on_mouse_message)

    def _mouse_move(self, event_data):
        if self._data is None:
            # Nothing to display, so exit
            return

        xc = event_data['domain']['x']
        yc = event_data['domain']['y']

        # get the array indices into the data so that we can get data values
        x_index = int(np.floor(xc + 0.5))
        y_index = int(np.floor(yc + 0.5))

        # Check that the index is in the array.
        in_image = (self._data.shape[1] > x_index >= 0) and (self._data.shape[0] > y_index >= 0)
        if in_image:
            val = self._data[y_index, x_index]
        else:
            val = None

        if val is not None:
            value = f'value: {val:8.3f}'
        else:
            value = 'value: N/A'

        pixel_location = f'X:  {xc:.2f}  Y:  {yc:.2f}'
        if self._wcs is not None:
            sky = self._wcs.pixel_to_world(yc, xc)
            ra_dec = f'RA: {sky.icrs.ra:3.7f} Dec: {sky.icrs.dec:3.7f}'
        else:
            ra_dec = ''
        self._cursor.value = ', '.join([pixel_location, ra_dec, value])

    def _mouse_click(self, event_data):
        if self._data is None:
            # Nothing to display, so exit
            return

        xc = event_data['domain']['x']
        yc = event_data['domain']['y']

        if self.click_center:
            self.center_on((xc, yc))

        if self.is_marking:
            print('marky marking')
            # Just hand off to the method that actually does the work
            self._add_new_single_marker(xc, yc)

    def _add_new_single_marker(self, x_mark, y_mark):
        # We have location of the new marker and should have the name
        # of the marker tag and the marker style, so just need to update
        # the table and draw the new maker.

        marker_name = self._marker_table._interactive_marker_set_name
        # update the marker table
        self._marker_table.add_markers([x_mark], [y_mark], marker_name=marker_name)

        # First approach: get any current markers by that name, add this one
        # remove the old ones and draw the new ones.
        marks = self.get_markers(marker_name=marker_name)
        self._astro_im.plot_named_markers(marks['x'], marks['y'],
                                          marker_name,
                                          color=self.marker['color'],
                                          size=self.marker['radius']**2,
                                          style=self.marker['type'])

    # Update the viewport to match changes in the UI
    def _init_watch_image_changes(self):
        """
        Watch for changes to the image scale, which indicate the user
        has either changed the zoom or has panned, and update the zoom_level.
        """
        def update_zoom_level(event):
            """
            Watch for changes in the zoom level from the viewer.
            """

            old_width = self.get_viewport(sky_or_pixel='pixel', image_label=self._current_image_label)["fov"]
            new_width = self._astro_im.get_current_width()
            if new_width is None or self._updating_viewport:
                # There is no image yet, or this object is in the process
                # of changing the zoom, so return
                return

            # Do nothing if the zoom has not changed
            if np.abs(new_width - old_width) > 1e-3:
                # Let the zoom_level handler know the GUI itself
                # generated this zoom change which means the GUI
                # does not need to be updated.
                self._viewport_change_source_is_gui = True
                self.set_viewport(fov=new_width)

        # Observe changes to the maximum of the x scale. Observing the y scale
        # or the minimum instead of the maximum is also fine.
        x_scale = self._astro_im._scales['x']

        # THIS IS TERRIBLE AND MAKES THINGS SUPER LAGGY!!!! Needs to be
        # throttled or something. Look at the ImageGL observe options.
        x_scale.observe(update_zoom_level, names='max')

    def _interval_and_stretch(self, stretch=None, cuts=None):
        """
        Stretch and normalize the data before sending to the viewer.
        """
        interval = cuts or self._default_cuts
        intervaled = interval(np.asarray(self._data))

        stretch = stretch or self._default_stretch
        if stretch:
            stretched = stretch(intervaled)
        else:
            stretched = intervaled

        return stretched

    def _send_data(self, reset_view=True, stretch=None, cuts=None):
        if self._data is not None:
            self._astro_im.set_data(self._interval_and_stretch(stretch=stretch, cuts=cuts),
                                    reset_view=reset_view)
        #self.zoom_level = self._astro_im.get_zoom_level()

    @property
    def _current_image_label(self):
        """
        Image label for the most recently loaded image
        """
        return list(self._images.keys())[-1]

    def set_stretch(self, value, image_label=None, **kwargs):
        super().set_stretch(value, image_label=image_label, **kwargs)
        self._send_data(stretch=value)

    def set_cuts(self, value, image_label=None, **kwargs):
        super().set_cuts(value, image_label=image_label, **kwargs)
        self._send_data(cuts=self.get_cuts(image_label=image_label))

    @property
    def viewer(self):
        return self._astro_im

    # The methods, grouped loosely by purpose
    def load_image(self, image, image_label=None, **kwargs):
        super().load_image(image, image_label=image_label, **kwargs)
        data = self.get_image(image_label=image_label)

        self._data = data.data if isinstance(data, NDData) else data
        self._send_data()

    # Saving contents of the view and accessing the view
    def save(self, filename, overwrite=False, **kwargs):
        """
        Saving for this backend requires a running browser.
        """
        p = Path(filename)

        if p.suffix == '.png':
            self._astro_im.save_png(filename, overwrite=overwrite)
        elif p.suffix == '.svg':
            self._astro_im.save_svg(filename, overwrite=overwrite)
        else:
            raise ValueError('Saving is not supported for that'
                             'file type. Use .png or .svg')

    def set_colormap(self, cmap_name, image_label=None, **kwargs):
        super().set_colormap(cmap_name, image_label=image_label, **kwargs)
        self._astro_im.set_color(bqcolors(cmap_name, reverse=False))

    def load_catalog(
        self,
        table,
        **kwargs
    ):
        super().load_catalog(table, **kwargs)
        catalog_label = kwargs.pop("catalog_label", None)
        this_catalog = self.get_catalog(catalog_label=catalog_label)
        self._astro_im.plot_named_markers(
            this_catalog["x"],
            this_catalog["y"],
            str(catalog_label),
            **self.get_catalog_style(catalog_label=catalog_label)
        )

    def set_catalog_style(
            self,
            catalog_label=None,
            shape="circle",
            color="red",
            size=5,
            **kwargs
        ):
        super().set_catalog_style(
            catalog_label=catalog_label,
            shape=shape,
            color=color,
            size=size,
            **kwargs
        )
        this_catalog = self.get_catalog(catalog_label=catalog_label)

        self._astro_im.plot_named_markers(
            this_catalog["x"],
            this_catalog["y"],
            str(catalog_label),
            color=color,
            size=size**2,  # bqplot expects size in pixels squared
            shape=shape,
        )

    def remove_catalog(self, catalog_label=None, **kwargs):
        super().remove_catalog(catalog_label, **kwargs)
        if catalog_label == "*":
            # Remove all catalogs
            self._astro_im.remove_markers()
        else:
            self._astro_im.remove_named_markers(str(catalog_label))

    def set_viewport(self, center=None, fov=None, image_label=None, **kwargs):
        # This will handle all of the WCS stuff, which we need in the event
        # the fov is in sky coordinates.
        super().set_viewport(center=center, fov=fov, image_label=image_label, **kwargs)

        # Get the viewport in pixel coordinates
        viewport = self.get_viewport(image_label=image_label, sky_or_pixel='pixel')
        if not self._viewport_change_source_is_gui:
            # This is used in the viewport handler to suppress
            # handling during a programmatic change to the
            # viewport.
            self._updating_viewport = True
            # Set the center of the image to the center of the viewport
            # Note the coordinates are reversed because that is how it goes
            # with image display vs array indices.
            self._astro_im.center = viewport['center']

            # Set the size of the image to the size of the viewport
            self._astro_im.set_size(viewport['fov'], direction='x')
            self._updating_viewport = False
        else:
            # The GUI is the source of the change, so do not update the
            # image center or size.
            self._viewport_change_source_is_gui = False

    def get_viewport(self, sky_or_pixel=None, image_label=None, **kwargs):
        return super().get_viewport(sky_or_pixel, image_label, **kwargs)

    @property
    def print_out(self):
        """
        Return an output widget for display in the notebook which
        captures any printed output produced by the viewer widget.

        Intended primarily for debugging.
        """
        return self._print_out
