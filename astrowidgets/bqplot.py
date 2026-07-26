from contextlib import ExitStack, contextmanager
from pathlib import Path

import numpy as np
from astropy.nddata import NDData
import astropy.visualization as apviz

from bqplot import Figure, LinearScale, Axis, ColorScale, PanZoom, Scatter
from bqplot_image_gl import ImageGL
from bqplot_image_gl.interacts import (MouseInteraction,
                                       keyboard_events, mouse_events)

import ipywidgets as ipw

from matplotlib import pyplot
from matplotlib.colors import to_hex

from astro_image_display_api.image_viewer_logic import (
    ImageViewerLogic,
    docs_from_image_viewer_logic_if_missing,
)

from astrowidgets.cursor_info import CursorInfoMixin


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

        scale_x = LinearScale(min=0, max=1,  # self._image_shape[1],
                              allow_padding=False)
        scale_y = LinearScale(min=0, max=1,  # self._image_shape[0],
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
        self.center = (self._image_shape[1] / 2,
                       self._image_shape[0] / 2)

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
                scale_factor = 1 / figure_aspect
            else:
                scale_factor = figure_aspect

            self._scales[reset_scale].min = 0
            self._scales[reset_scale].max = frozen_width[reset_scale] * scale_factor
            self.center = current_center

    @contextmanager
    def _hold_all_sync(self):
        """
        Batch trait sync messages for the image mark, its color scale and
        both axis scales so that the front end receives a single state
        update per widget, and hence redraws once, instead of redrawing
        after every trait assignment.

        These are separate widgets, so each still syncs its own message and
        the front end redraws on each. On exit the contexts release -- and
        each widget flushes -- in reverse of the order entered here, so the
        entry order is chosen to make the image mark flush first, then the
        color scale, then the scales. That way the front end never draws
        the old image against the new scales (a "refit" flash) nor recolors
        the old image before the new data arrives.
        """
        with self._scales['y'].hold_sync(), self._scales['x'].hold_sync(), \
                self._image.scales['image'].hold_sync(), \
                self._image.hold_sync():
            yield

    def set_data(self, image_data, reset_view=True):
        self._image_shape = image_data.shape

        with self._hold_all_sync():
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
            scale_to_set.min = cen[direction] - size / 2
            scale_to_set.max = cen[direction] + size / 2

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
        sc = Scatter(scales=scale_dict,
                     x=np.asarray(x), y=np.asarray(y),
                     colors=[color],
                     default_size=size,
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
@docs_from_image_viewer_logic_if_missing
class ImageWidget(ipw.VBox, CursorInfoMixin, ImageViewerLogic):
    def __init__(self, *args, display_width=500, display_aspect_ratio=1):
        super().__init__(*args)
        # ImageViewerLogic is a dataclass; we do not run its __init__, so run
        # the post-init hook it would otherwise provide to set up its state.
        ImageViewerLogic.__post_init__(self)

        self._astro_im = _AstroImage(display_width=display_width,
                                     viewer_aspect_ratio=display_aspect_ratio)
        # Cut out the sky background at the bottom and clip only the
        # brightest pixels at the top.
        self._default_cuts = apviz.AsymmetricPercentileInterval(30, 96)
        self._default_stretch = None
        self._default_colormap = 'Greys_r'
        # The API layer stores settings per image, so nothing can be stored
        # before a load; apply the default colormap directly to the front
        # end so the empty viewer already shows it. When the first image is
        # loaded, _render_image stores the widget defaults for it so that
        # the get_* methods report what is displayed.
        self._astro_im.set_color(bqcolors(self._default_colormap))

        self._data = None
        self._wcs = None

        # Use this to manage whether or not to send changes in zoom level
        # to the viewer.
        self._viewport_change_source_is_gui = False

        # Guards re-entrancy while we programmatically update the viewport.
        self._updating_viewport = False

        # While True, _refresh_display does nothing; set by _defer_refresh.
        self._refresh_deferred = False

        # Provide an Output widget to which prints can be directed for
        # debugging.
        self._print_out = ipw.Output()

        self._init_mouse_callbacks()
        self._init_watch_image_changes()
        self.children = [self._astro_im, self._init_cursor_info()]

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
        self._update_cursor_text(event_data['domain']['x'],
                                 event_data['domain']['y'])

    def _mouse_click(self, event_data):
        if self._data is None:
            # Nothing to display, so exit
            return

        # Interactive click features (click-to-center, interactive marking)
        # were removed in the migration to the astro-image-display-api and
        # will be rebuilt, see #201. Until then this is intentionally a
        # no-op so that clicks do not raise and user-registered on_msg
        # callbacks still run.

    # Update the viewport to match changes in the UI
    def _init_watch_image_changes(self):
        """
        Watch for changes to the image scales, which indicate the user
        has either changed the zoom or has panned, and update the stored
        viewport (center and field of view) to match.
        """
        def update_viewport(event):
            """
            Watch for changes in the viewport (pan or zoom) from the viewer.
            """
            new_width = self._astro_im.get_current_width()
            if (new_width is None or self._updating_viewport
                    or not self._displayed_image_labels):
                # There is no image yet, or this object is in the process
                # of changing the viewport, so return
                return

            # A pan or zoom in the GUI is always a change to the view of
            # the displayed image.
            image_label = self._displayed_image_labels[0]
            current = self.get_viewport(sky_or_pixel='pixel',
                                        image_label=image_label)
            old_width = current["fov"]
            old_center = current["center"]
            new_center = self._astro_im.center

            width_changed = np.abs(new_width - old_width) > 1e-3
            center_changed = (
                old_center is None
                or np.abs(new_center[0] - old_center[0]) > 1e-3
                or np.abs(new_center[1] - old_center[1]) > 1e-3
            )

            # Do nothing if neither the zoom nor the pan has changed
            if width_changed or center_changed:
                # Let the viewport handler know the GUI itself generated this
                # change, which means the GUI does not need to be updated.
                self._viewport_change_source_is_gui = True
                self.set_viewport(center=new_center, fov=new_width,
                                  image_label=image_label)

        # Observe changes to the maximum of both the x and y scales so that
        # horizontal pan, vertical pan, and zoom are all detected. Observing
        # the maximum (rather than the minimum) of each scale is sufficient
        # because a pan shifts both ends of a scale.
        #
        # If things seem laggy in the future, check whether throttling
        # the updates helps.
        for scale in (self._astro_im._scales['x'], self._astro_im._scales['y']):
            scale.observe(update_viewport, names='max')

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

    def _refresh_display(self, image_label=None):
        """
        Recompute the displayed array from the cuts and stretch stored
        for the image and send it to the viewer. The viewport (zoom/pan)
        is left untouched; it is owned by the ``_apply_viewport`` hook.
        """
        if self._data is None or self._refresh_deferred:
            return

        self._send_data(cuts=self.get_cuts(image_label=image_label),
                        stretch=self.get_stretch(image_label=image_label),
                        reset_view=False)

    @contextmanager
    def _defer_refresh(self):
        """
        Make _refresh_display a no-op inside the block so that several
        settings can be stored with a single recomputation of the
        displayed array. The caller refreshes after the block.
        """
        self._refresh_deferred = True
        try:
            yield
        finally:
            self._refresh_deferred = False

    # ------------------------------------------------------------------
    # Rendering hooks
    #
    # The public API methods inherited from ImageViewerLogic are templates
    # that own all state handling and label resolution, then call these
    # hooks with already-resolved labels to push the stored state into the
    # bqplot figure. The _apply_* hooks are only called for the displayed
    # image.
    # ------------------------------------------------------------------
    @contextmanager
    def _batch_update(self):
        """
        Batch the front-end updates of a group of state changes.

        Overrides the no-op hook from
        `~astro_image_display_api.ImageViewerLogic`; ``load_image`` wraps
        its state changes and rendering-hook calls in this context.

        The composition is: hold the front-end sync of the image mark and
        the scales for the whole block (so each widget sends one state
        message, with the image mark flushing first -- see
        ``_AstroImage._hold_all_sync``), and defer the display refreshes
        that the ``_apply_*`` hooks trigger, collapsing them into a single
        recomputation of the displayed array that runs while the front-end
        sync is still held.
        """
        with ExitStack() as stack:
            # Entered first, so it releases last: everything below happens
            # while the front-end sync is held.
            stack.enter_context(self._astro_im._hold_all_sync())
            # Runs on exit, after _defer_refresh has released: the one
            # refresh that replaces the refreshes deferred in the block.
            stack.callback(self._refresh_after_batch)
            stack.enter_context(self._defer_refresh())
            yield

    def _refresh_after_batch(self):
        """
        Recompute and send the displayed array once after a batched update.

        The viewport is not touched here: the ``_apply_viewport`` hook has
        already pushed the stored viewport into the scales during the batch.
        """
        for image_label in self._displayed_image_labels:
            self._refresh_display(image_label=image_label)

    def _render_image(self, image_label):
        """
        Make the image stored under a label the one the widget displays.

        Overrides the no-op rendering hook from
        `~astro_image_display_api.ImageViewerLogic`. ``load_image`` calls it
        after storing the image data and settings; the ``_apply_*`` hooks
        are called afterwards to push the stored settings into the display.

        Parameters
        ----------
        image_label : str
            Resolved label of the image to display.
        """
        info = self._images[image_label]

        if info.colormap is None:
            # load_image carries the displayed image's cuts, stretch and
            # colormap forward to the new image, and any image that has been
            # displayed has a colormap (this very block guarantees it), so a
            # missing colormap means nothing was carried forward: this image
            # is the first to be displayed. Store the widget's defaults for
            # it so that the get_* methods report what is displayed.
            info.colormap = self._default_colormap
            info.cuts = self._default_cuts
            if self._default_stretch is not None:
                info.stretch = self._default_stretch

        data = info.data
        self._data = data.data if isinstance(data, NDData) else data
        self._wcs = info.wcs

    def _apply_cuts(self, image_label):
        """
        Re-display the image using the cut levels stored for a label.

        Overrides the no-op rendering hook from
        `~astro_image_display_api.ImageViewerLogic`; only called when the
        label is displayed. Changing the cuts only affects the color
        mapping, so the current viewport (zoom/pan) is left untouched.

        Parameters
        ----------
        image_label : str
            Resolved label whose stored cuts to apply.
        """
        self._refresh_display(image_label=image_label)

    def _apply_stretch(self, image_label):
        """
        Re-display the image using the stretch stored for a label.

        Overrides the no-op rendering hook from
        `~astro_image_display_api.ImageViewerLogic`; only called when the
        label is displayed. Changing the stretch only affects the color
        mapping, so the current viewport (zoom/pan) is left untouched.

        Parameters
        ----------
        image_label : str
            Resolved label whose stored stretch to apply.
        """
        self._refresh_display(image_label=image_label)

    def _apply_colormap(self, image_label):
        """
        Push the colormap stored for a label into the bqplot color scale.

        Overrides the no-op rendering hook from
        `~astro_image_display_api.ImageViewerLogic`; only called when the
        label is displayed.

        Parameters
        ----------
        image_label : str
            Resolved label whose stored colormap to apply.
        """
        self._astro_im.set_color(
            bqcolors(self.get_colormap(image_label=image_label)))

    def _apply_viewport(self, image_label):
        """
        Push the viewport stored for a label into the bqplot scales.

        Overrides the no-op rendering hook from
        `~astro_image_display_api.ImageViewerLogic`; only called when the
        label is displayed.

        Parameters
        ----------
        image_label : str
            Resolved label whose stored viewport (center and field of view)
            to apply.
        """
        if self._viewport_change_source_is_gui:
            # The GUI itself (a pan or zoom in the browser) generated this
            # viewport change, so the front end is already up to date.
            self._viewport_change_source_is_gui = False
            return

        # Get the viewport in pixel coordinates; the API layer handles all
        # of the WCS stuff in the event the stored fov is in sky units.
        viewport = self.get_viewport(image_label=image_label,
                                     sky_or_pixel='pixel')

        # Suppress the handling of the scale changes made below as if they
        # were a pan or zoom from the GUI.
        self._updating_viewport = True
        try:
            self._astro_im.center = viewport['center']
            self._astro_im.set_size(viewport['fov'], direction='x')
        finally:
            self._updating_viewport = False

    def _draw_catalog(self, catalog_label):
        """
        Draw (or redraw) a catalog's markers as a bqplot scatter mark.

        Overrides the no-op rendering hook from
        `~astro_image_display_api.ImageViewerLogic`; called by
        ``load_catalog`` and ``set_catalog_style``.

        Parameters
        ----------
        catalog_label : str
            Resolved label of the catalog to draw. Its markers are drawn
            under a mark id equal to the label, replacing any previous
            scatter mark for that catalog.
        """
        catalog = self.get_catalog(catalog_label=catalog_label)
        style = self.get_catalog_style(catalog_label=catalog_label)

        self._astro_im.plot_named_markers(
            catalog["x"],
            catalog["y"],
            catalog_label,
            color=style.get("color", "red"),
            # bqplot expects the size in pixels squared
            size=style.get("size", 5)**2,
            shape=style.get("shape", "circle"),
        )

    def _remove_catalog_marks(self, catalog_label):
        """
        Remove a catalog's scatter mark from the bqplot figure.

        Overrides the no-op rendering hook from
        `~astro_image_display_api.ImageViewerLogic`; ``remove_catalog``
        calls it once per removed catalog (after expanding ``"*"``).

        Parameters
        ----------
        catalog_label : str
            Resolved label of the catalog whose markers to remove.
        """
        self._astro_im.remove_named_markers(catalog_label)

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

    @property
    def print_out(self):
        """
        Return an output widget for display in the notebook which
        captures any printed output produced by the viewer widget.

        Intended primarily for debugging.
        """
        return self._print_out
