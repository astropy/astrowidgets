from pathlib import Path

import numpy as np

from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.io import fits
from astropy.nddata import CCDData, NDData
from astropy.table import Table
from astropy import units as u
import astropy.visualization as apviz

from bqplot import Figure, LinearScale, Axis, ColorScale, PanZoom, ScatterGL
from bqplot_image_gl import ImageGL
from bqplot_image_gl.interacts import (MouseInteraction,
                                       keyboard_events, mouse_events)

import ipywidgets as ipw

from matplotlib import pyplot
from matplotlib.colors import to_hex

from astro_image_display_api import ImageViewerInterface

import traitlets as trait


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

        self._image = ImageGL(scales=scales_image)

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
    def center(self):
        """
        Center of current view in pixels in x, y.
        """
        x_center = (self._scales['x'].min + self._scales['x'].max) / 2
        y_center = (self._scales['y'].min + self._scales['y'].max) / 2
        return (x_center, y_center)

    @property
    def scale_widths(self):
        width_x = self._scales['x'].max - self._scales['x'].min
        width_y = self._scales['y'].max - self._scales['y'].min
        return (width_x, width_y)

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
        scale_to_set = self._scales[direction]
        cen = {}
        cen['x'], cen['y'] = self.center
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

    def get_zoom_level(self):
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
        figure_width = float(self._figure.layout.width[:-2])
        scale_width = self.scale_widths[1]

        return figure_width / scale_width

    def plot_named_markers(self, x, y, mark_id, color='yellow',
                           size=100, style='circle'):
        scale_dict = dict(x=self._scales['x'], y=self._scales['y'])
        sc = ScatterGL(scales=scale_dict,
                       x=x, y=y,
                       colors=[color],
                       default_size=100,
                       marker=style,
                       fill=False)

        self._scatter_marks[mark_id] = sc
        self._update_marks()

    def remove_named_markers(self, mark_id):
        if isinstance(mark_id, str):
            if mark_id == 'all':
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


class MarkerTableManager:
    """
    Table for keeping track of positions and names of sets of
    logically-related markers.
    """
    def __init__(self):
        # These column names are for internal use.
        self._xcol = 'x'
        self._ycol = 'y'
        self._names = 'marker name'
        self._marktags = set()
        # Let's have a default name for the tag too:
        self.default_mark_tag_name = 'default-marker-name'
        self._interactive_marker_set_name_default = 'interactive-markers'
        self._interactive_marker_set_name = self._interactive_marker_set_name_default
        self._init_table()

    def _init_table(self):
        self._table = Table(names=(self._xcol, self._ycol, self._names),
                            dtype=(np.float64, np.float64, 'str'))

    @property
    def xcol(self):
        return self._xcol

    @property
    def ycol(self):
        return self._ycol

    @property
    def names(self):
        return self._names

    @property
    def marker_names(self):
        return sorted(set(self._table[self.names]))

    def add_markers(self, x_mark, y_mark,
                    marker_name=None):

        if marker_name is None:
            marker_name = self.default_mark_tag_name

        self._marktags.add(marker_name)
        for x, y in zip(x_mark, y_mark):
            self._table.add_row([x, y, marker_name])

    def get_markers(self, marker_name):
        if isinstance(marker_name, str):
            if marker_name == 'all':
                return self._table.copy()
            else:
                marker_name = [marker_name]

        row_mask = np.zeros(len(self._table), dtype=bool)
        for name in marker_name:
            row_mask |= self._table[self._names] == name

        return self._table[row_mask]

    def get_all_markers(self):
        return self._table.copy()

    def remove_markers(self, marker_name):
        if isinstance(marker_name, str):
            if marker_name == 'all':
                self.remove_all_markers()
                return
            else:
                marker_name = [marker_name]

        to_delete = np.zeros(len(self._table), dtype=bool)
        for name in marker_name:
            to_delete |= self._table[self._names] == name

        # Only keep the things that don't match
        self._table = self._table[~to_delete]

    def remove_all_markers(self):
        self._init_table()


"""
next(iter(imviz.app._viewer_store.values())).figure
"""
STRETCHES = dict(
    linear=apviz.LinearStretch,
    sqrt=apviz.SqrtStretch,
    histeq=apviz.HistEqStretch,
    log=apviz.LogStretch
    # ...
)


class ImageWidget(ipw.VBox):
    click_center = trait.Bool(default_value=False).tag(sync=True)
    click_drag = trait.Bool(default_value=False).tag(sync=True)
    scroll_pan = trait.Bool(default_value=False).tag(sync=True)
    image_width = trait.Int(help="Width of the image (not viewer)").tag(sync=True)
    image_height = trait.Int(help="Height of the image (not viewer)").tag(sync=True)
    zoom_level = trait.Float(help="Current zoom of the view").tag(sync=True)
    marker = trait.Any(help="Markers").tag(sync=True)

    # Wonder if there is a way to force this to be a subclass of BaseInterval?
    # Yes, with Instance (or maybe Type?)
    cuts = trait.Any(help="Cut levels", allow_none=True).tag(sync=False)
    cursor = "foo"

    # Wonder if there is a way to force this to be a subclass of BaseStretch?
    # Yes, with Instance (or maybe Type?)
    stretch = trait.Unicode(help='Stretch algorithm name', allow_none=True).tag(sync=True)
    stretch_options = trait.List(trait.Unicode(),
                                 default_value=list(STRETCHES.keys()),
                                 help='Stretch algorithm options').tag(sync=True)

    autocut_options = trait.List(trait.Unicode(),
                                 default_value=['minmax', 'percentile'],
                                 help='Autocut options').tag(sync=True)

    def __init__(self, *args, image_width=500, image_height=500):
        super().__init__(*args)
        self.RESERVED_MARKER_SET_NAMES = ['all']
        self.image_width = image_width
        self.image_height = image_height
        viewer_aspect = self.image_width / self.image_height
        self._astro_im = _AstroImage(display_width=self.image_width,
                                     viewer_aspect_ratio=viewer_aspect)
        self._interval = None
        self._stretch = None
        self.set_colormap('Greys_r')
        self._marker_table = MarkerTableManager()
        self._data = None
        self._wcs = None

        # Marker state is going away
        self._is_marking = False

        self.scroll_pan = True
        # Use this to manage whether or not to send changes in zoom level
        # to the viewer.
        self._zoom_source_is_gui = False
        # Use this to let the method monitoring changes coming from the
        # image know that the ImageWidget itself is in the process of
        # updating the zoom.
        self._updating_zoom = False

        # Provide an Output widget to which prints can be directed for
        # debugging.
        self._print_out = ipw.Output()

        self.marker = {'color': 'red', 'radius': 20, 'type': 'square'}
        self.cuts = apviz.AsymmetricPercentileInterval(1, 99)

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

    ## THIS CAN BE DROPPED, I THINK -- REPLACE WITH VIEWPORT
    def _init_watch_image_changes(self):
        """
        Watch for changes to the image scale, which indicate the user
        has either changed the zoom or has panned, and update the zoom_level.
        """
        def update_zoom_level(event):
            """
            Watch for changes in the zoom level from the viewer.
            """

            old_zoom = self.zoom_level
            new_zoom = self._astro_im.get_zoom_level()
            if new_zoom is None or self._updating_zoom:
                # There is no image yet, or this object is in the process
                # of changing the zoom, so return
                return

            # Do nothing if the zoom has not changed
            if np.abs(new_zoom - old_zoom) > 1e-3:
                # Let the zoom_level handler know the GUI itself
                # generated this zoom change which means the GUI
                # does not need to be updated.
                self._zoom_source_is_gui = True
                self.zoom_level = new_zoom

        # Observe changes to the maximum of the x scale. Observing the y scale
        # or the minimum instead of the maximum is also fine.
        x_scale = self._astro_im._scales['x']

        # THIS IS TERRIBLE AND MAKES THINGS SUPER LAGGY!!!! Needs to be
        # throttled or something. Look at the ImageGL observe options.
        x_scale.observe(update_zoom_level, names='max')

    def _interval_and_stretch(self):
        """
        Stretch and normalize the data before sending to the viewer.
        """
        interval = self._get_interval()
        intervaled = interval(self._data)

        stretch = self._get_stretch()
        if stretch:
            stretched = stretch(intervaled)
        else:
            stretched = intervaled

        return stretched

    def _send_data(self, reset_view=True):
        self._astro_im.set_data(self._interval_and_stretch(),
                                reset_view=reset_view)
        self.zoom_level = self._astro_im.get_zoom_level()

    def _get_interval(self):
        if self._interval is None:
            return apviz.MinMaxInterval()
        else:
            return self._interval

    def _get_stretch(self):
        return self._stretch

    @trait.validate('stretch')
    def _validate_stretch(self, proposal):
        proposed_stretch = proposal['value']
        if (proposed_stretch not in STRETCHES.keys() and
            proposed_stretch is not None):

            raise ValueError(f'{proposed_stretch} is not a valid value. '
                              'The stretch must be one of None or '
                              'one of these values: '
                             f'{sorted(STRETCHES.keys())}')

        return proposed_stretch

    @trait.observe('stretch')
    def _observe_stretch(self, change):
        if change['new'] == 'histeq':
            self._stretch = STRETCHES[change['new']](self._data)
        else:
            self._stretch = STRETCHES[change['new']]() if change['new'] else None

        if self._stretch is not None and change['new'] != change['old']:
            self._send_data()

    @trait.validate('cuts')
    def _validate_cuts(self, proposal):
        # Allow these:
        # - a two-item thing (tuple, list, whatever)
        # - an Astropy interval
        # - a name, like 'minmax'
        # - None
        proposed_cuts = proposal['value']

        bad_value_error = (f"{proposed_cuts} is not a valid value. "
                           "cuts must be one of None, "
                           "an astropy interval, name, or list/tuple "
                           "of length 2.")

        if ((proposed_cuts is None) or
            isinstance(proposed_cuts, apviz.BaseInterval)):
            return proposed_cuts
        elif isinstance(proposed_cuts, str):
            if proposed_cuts not in self.autocut_options:
                raise ValueError(bad_value_error)

            match proposed_cuts:
                case 'minmax':
                    proposed_cuts = apviz.MinMaxInterval().get_limits(self._data)
                case 'percentile':
                    proposed_cuts = apviz.PercentileInterval(1, 99)
                case _:
                    raise ValueError(bad_value_error)
            return proposed_cuts
        else:
            try:
                length = len(proposed_cuts)
                if length != 2:
                    raise ValueError('Cut levels must be given as (low, high).'
                                        + bad_value_error)

                # Tests expect this to be a tuple...
                proposed_cuts = tuple(proposed_cuts)
            except (TypeError, AssertionError):
                raise ValueError(bad_value_error)

            return proposed_cuts

    @trait.observe('cuts')
    def _observe_cuts(self, change):
        # This needs to handle only the case when the cuts is a
        # tuple/list of length 2. That is interpreted as a ManualInterval.
        cuts = change['new']
        if cuts is not None:
            if not isinstance(cuts, apviz.BaseInterval):
                self._interval = apviz.ManualInterval(*cuts)
            else:
                self._interval = cuts
        if self._data is not None:
            self._send_data(reset_view=False)

    @trait.observe('zoom_level')
    def _update_zoom_level(self, change):
        zl = change['new']
        if not self._zoom_source_is_gui:
            # User has changed the zoom value so update the viewer
            self._updating_zoom = True
            self._astro_im.set_zoom_level(zl)
            self._updating_zoom = False
        else:
            # GUI updated the value so do nothing except reset the source
            # of the event
            self._zoom_source_is_gui = False

    ## ✂️ DROP THIS FROM HERE ✂️ to next scissors
    def _currently_marking_error_msg(self, caller):
        return (f'Cannot set {caller} while doing interactive '
                f'marking. Call the stop_marking() method to '
                f'stop marking and then set {caller}.')

    @trait.validate('click_drag')
    def _validate_click_drag(self, proposal):
        cd = proposal['value']
        if cd and self._is_marking:
            raise ValueError(self._currently_marking_error_msg('click_drag'))
        return cd

    @trait.observe('click_drag')
    def _update_viewer_pan(self, change):
        # Turn of click-to-center
        if change['new']:
            self.click_center = False

        self._astro_im.set_pan(change['new'])

    @trait.observe('scroll_pan')
    def _update_viewer_zoom_scroll(self, change):
        self._astro_im.set_scroll_zoom(change['new'])

    @trait.validate('click_center')
    def _validate_click_center(self, proposal):
        new = proposal['value']
        if new and self._is_marking:
            raise ValueError(self._currently_marking_error_msg('click_center'))
        return new

    @trait.observe('click_center')
    def _update_click_center(self, change):
        if change['new']:
            # click_center has been turned on, so turn off click_drag
            self.click_drag = False

    @trait.observe('cursor')
    def _update_cursor_position(self, change):
        if change['new'] == 'top':
            self.layout.flex_flow = 'column-reverse'
            self._cursor.layout.visibility = 'visible'
        elif change['new'] == 'bottom':
            self.layout.flex_flow = 'column'
            self._cursor.layout.visibility = 'visible'
        elif change['new'] is None:
            self._cursor.layout.visibility = 'hidden'

    # ✂️ END THIS CUT

    @property
    def viewer(self):
        return self._astro_im

    # CUT
    @property
    def is_marking(self):
        """`True` if in marking mode, `False` otherwise.
        Marking mode means a mouse click adds a new marker.
        This does not affect :meth:`add_markers`.

        """
        return self._is_marking

    # CUT
    @property
    def _default_mark_tag_name(self):
        """
        This is only here to make a test pass -- it should probably either
        be part of the API or not tested.
        """
        return self._marker_table.default_mark_tag_name

    # The methods, grouped loosely by purpose

    # Methods for loading data
    def load(self, data):
        if isinstance(data, NDData):
            self._load_nddata(data)
        elif isinstance(data, (str, Path, fits.ImageHDU, fits.CompImageHDU, fits.PrimaryHDU)):
            self._load_fits(data)
        else:
            self._load_array(data)

    def _load_fits(self, file_name_or_HDU, reset_view=True):
        if isinstance(file_name_or_HDU, (str, Path)):
            ccd = CCDData.read(file_name_or_HDU)
        elif isinstance(file_name_or_HDU,
                        (fits.ImageHDU, fits.CompImageHDU, fits.PrimaryHDU)):
            try:
                ccd_unit = u.Unit(file_name_or_HDU.header['bunit'])
            except (KeyError, ValueError):
                ccd_unit = u.dimensionless_unscaled
            ccd = CCDData(file_name_or_HDU.data,
                          header=file_name_or_HDU.header,
                          unit=ccd_unit)
        else:
            raise ValueError(f'{file_name_or_HDU} is an invalid value. It must'
                             ' be a string or an astropy.io.fits HDU.')

        self._ccd = ccd
        self._data = ccd.data
        self._wcs = ccd.wcs
        self._send_data(reset_view=reset_view)

    def _load_array(self, array, reset_view=True):
        self._data = array
        self._send_data(reset_view=reset_view)

    def _load_nddata(self, data, reset_view=True):
        self._ccd = data
        self._data = self._ccd.data
        self._wcs = data.wcs

        self._send_data(reset_view=reset_view)

    # Saving contents of the view and accessing the view
    def save(self, filename, overwrite=False):
        p = Path(filename)

        if p.suffix == '.png':
            self._astro_im.save_png(filename, overwrite=overwrite)
        elif p.suffix == '.svg':
            self._astro_im.save_svg(filename, overwrite=overwrite)
        else:
            raise ValueError('Saving is not supported for that'
                             'file type. Use .png or .svg')

    def set_colormap(self, cmap_name, reverse=False):
        self._astro_im.set_color(bqcolors(cmap_name, reverse=reverse))
        self._colormap = cmap_name

    @property
    def colormap_options(self):
        return pyplot.colormaps()

    def _validate_marker_name(self, marker_name):
        if marker_name in self.RESERVED_MARKER_SET_NAMES:
            raise ValueError(
                f"The marker name {marker_name} is not allowed. Any name is "
                f"allowed except these: {', '.join(self.RESERVED_MARKER_SET_NAMES)}")

    def add_markers(self, table, x_colname='x', y_colname='y',
                    skycoord_colname='coord', use_skycoord=False,
                    marker_name=None):

        # Handle the case where marker_name is None
        if marker_name is None:
            marker_name = self._marker_table.default_mark_tag_name

        self._validate_marker_name(marker_name)

        if use_skycoord:
            if self._wcs is None:
                raise ValueError('The WCS for the image must be set to use '
                                 'world coordinates for markers.')

            x, y = self._wcs.world_to_pixel(table[skycoord_colname])
        else:
            x = table[x_colname]
            y = table[y_colname]

        # Update the table of marker names and positions
        self._marker_table.add_markers(x, y, marker_name=marker_name)

        # Update the figure itself, which expects all markers of
        # the same name to be plotted at once.
        marks = self.get_markers(marker_name)

        if marks:
            self._astro_im.plot_named_markers(marks['x'], marks['y'],
                                              marker_name,
                                              color=self.marker['color'],
                                              size=self.marker['radius']**2,
                                              style=self.marker['type'])

    def remove_markers(self, marker_name):
        # Remove from our tracking table
        self._marker_table.remove_markers(marker_name)

        # Remove from the visible canvas
        self._astro_im.remove_named_markers(marker_name)

    def reset_markers(self):
        self._marker_table.remove_all_markers()
        self._astro_im.remove_markers()

    def _prepare_return_marker_table(self, marks, x_colname='x', y_colname='y',
                                     skycoord_colname='coord'):
        if len(marks) == 0:
            return Table(names=['x', 'y', 'coord', 'marker name'])

        if (self._data is None) or (self._wcs is None):
            # Do not include SkyCoord column
            include_skycoord = False
        else:
            include_skycoord = True
            radec_col = []

        if include_skycoord:
            coords = self._wcs.pixel_to_world(marks[self._marker_table.xcol],
                                              marks[self._marker_table.ycol])
            marks[skycoord_colname] = coords

        # This might be a null op but should be harmless in that case
        marks.rename_column(self._marker_table.xcol, x_colname)
        marks.rename_column(self._marker_table.ycol, y_colname)

        return marks

    def get_marker_names(self):
        return self._marker_table.marker_names

    def get_markers(self, marker_name=None, x_colname='x', y_colname='y',
                            skycoord_colname='coord'):
        # We should always allow the default name. The case
        # where that table is empty will be handled in a moment.
        if isinstance(marker_name, str):
            if (marker_name not in self._marker_table.marker_names
                    and marker_name != self._marker_table.default_mark_tag_name
                    and marker_name != 'all'):
                return Table(names=["x", "y", "coord", "marker name"])

        if marker_name == "all":
            # This is a special case that is used in the tests
            # to get all markers.
            marks = self._marker_table.get_all_markers()
        else:
            if marker_name is None:
                # This is the default name
                marker_name = self.DEFAULT_MARKER_NAME
            marks = self._marker_table.get_markers(marker_name=marker_name)

        marks = self._prepare_return_marker_table(marks,
                                                  x_colname=x_colname,
                                                  y_colname=y_colname,
                                                  skycoord_colname=skycoord_colname)
        return marks

    def get_all_markers(self, x_colname='x', y_colname='y',
                        skycoord_colname='coord'):
        marks = self._marker_table.get_all_markers()
        marks = self._prepare_return_marker_table(marks,
                                                  x_colname=x_colname,
                                                  y_colname=y_colname,
                                                  skycoord_colname=skycoord_colname)
        return marks

    # Methods that modify the view
    def center_on(self, point):
        if isinstance(point, SkyCoord):
            if self._wcs is None:
                raise ValueError('The image must have a WCS to be able '
                                 'to center on a coordinate.')
            pixel = self._wcs.world_to_pixel(point)
        else:
            pixel = point

        self._astro_im.center = pixel

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
            raise u.UnitConversionError(f'dx is of type {dx_coord} but dy is of type {dy_coord} and they are not convertible.')


        if dx_coord == 'data':
            x, y = self._astro_im.center
            self.center_on((x + dx_val, y + dy_val))
        else:
            center_coord = self._wcs.pixel_to_world(*self._astro_im.center)
            # dx and dy in this case have units and we need to pass those units
            # in to offset.

            offset = SkyOffsetFrame(dx, dy, origin=center_coord.frame)
            new_center = SkyCoord(offset.transform_to(center_coord))

            # This is so much better only available in 4.3 or higher:
            # new_center = center_coord.spherical_offsets_by(dx_val, dy_val)

            self.center_on(new_center)

    def zoom(self, value):
        self.zoom_level = self.zoom_level * value

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
        mt = self._marker_table
        if marker_name is not None:
            self._validate_marker_name(marker_name)
            mt._interactive_marker_set_name = marker_name
        else:
            mt._interactive_marker_set_name = mt._interactive_marker_set_name_default
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
                self.reset_markers()

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
    def print_out(self):
        """
        Return an output widget for display in the notebook which
        captures any printed output produced by the viewer widget.

        Intended primarily for debugging.
        """
        return self._print_out
