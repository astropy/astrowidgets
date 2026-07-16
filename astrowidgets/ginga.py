"""AIDA-compliant image widget for Jupyter notebooks using the Ginga viewer."""

# STDLIB
import functools
from pathlib import Path
import warnings

# THIRD-PARTY
import numpy as np
from astropy.nddata import NDData
from astropy.table import Table
from astropy.visualization import (AsinhStretch, LinearStretch, LogStretch,
                                   PowerDistStretch, SinhStretch, SqrtStretch,
                                   SquaredStretch)

# Jupyter widgets
import ipywidgets as ipyw

# Ginga
from ginga import ColorDist
from ginga import cmap as ginga_cmap
from ginga.AstroImage import AstroImage
from ginga.canvas.CanvasObject import drawCatalog
from ginga.web.jupyterw.ImageViewJpw import EnhancedCanvasView
from ginga.util.wcs import ra_deg_to_str, dec_deg_to_str

from astro_image_display_api.image_viewer_logic import ImageViewerLogic

__all__ = ['ImageWidget']


class _AstropyStretchDist(ColorDist.ColorDistBase):
    """
    A ginga color distribution backed directly by an astropy stretch.

    Ginga has no native color-distribution class for some astropy stretches
    (e.g. ``PowerStretch`` (x**a), ``ContrastBiasStretch``, ``HistEqStretch``).
    Because every astropy stretch maps the unit interval onto itself, the
    stretch evaluated on ginga's 0..1 ramp *is* the lookup table ginga needs,
    so this adapter reproduces any astropy stretch exactly.
    """

    def __init__(self, hashsize, stretch, colorlen=None):
        self.stretch = stretch
        super().__init__(hashsize, colorlen=colorlen)

    def calc_hash(self):
        base = np.arange(0.0, float(self.hashsize), 1.0) / self.hashsize
        out = np.asarray(self.stretch(base, clip=True), dtype=float)
        out = np.clip(out, 0.0, 1.0)
        # normalize to color range
        ll = out * (self.colorlen - 1)
        self.hash = ll.astype(np.uint, copy=False)
        self.check_hash()

    def get_dist_pct(self, pct):
        # Inverse mapping, used to place color-bar ticks. Astropy stretches
        # expose ``.inverse``; fall back to the identity if one is unavailable.
        pct = np.asarray(pct, dtype=float)
        try:
            val = np.asarray(self.stretch.inverse(pct), dtype=float)
        except (NotImplementedError, AttributeError):
            val = pct
        return np.clip(val, 0.0, 1.0)

    def __str__(self):
        return type(self.stretch).__name__


def _ginga_dist_for_stretch(stretch, hashsize):
    """
    Build a ginga `~ginga.ColorDist.ColorDistBase` reproducing an astropy
    stretch, honoring its shape parameter.

    Ginga's parametrized distributions are reparametrizations of astropy's
    stretches, so where ginga has the family we configure its native class:
    the hyperbolic families match exactly, while log/power match the curve
    shape to within a quantization level (ginga normalizes by ``log(a)``/``a``
    where astropy uses ``log(a + 1)``/``a - 1``). Stretches ginga lacks --
    including ``PowerStretch`` (x**a), which is a different family than ginga's
    ``'power'`` (astropy's ``PowerDistStretch``) -- are handled exactly by
    `_AstropyStretchDist`.
    """
    if isinstance(stretch, AsinhStretch):
        factor = 1.0 / stretch.a
        return ColorDist.AsinhDist(hashsize, factor=factor,
                                   nonlinearity=np.arcsinh(factor))
    if isinstance(stretch, SinhStretch):
        factor = 1.0 / stretch.a
        return ColorDist.SinhDist(hashsize, factor=factor,
                                  nonlinearity=np.sinh(factor))
    if isinstance(stretch, LogStretch):
        return ColorDist.LogDist(hashsize, exp=stretch.a)
    if isinstance(stretch, PowerDistStretch):
        return ColorDist.PowerDist(hashsize, exp=stretch.a)
    # SquaredStretch subclasses PowerStretch, so check it before any plain
    # PowerStretch would fall through to the adapter.
    if isinstance(stretch, SquaredStretch):
        return ColorDist.SquaredDist(hashsize)
    if isinstance(stretch, SqrtStretch):
        return ColorDist.SqrtDist(hashsize)
    if (isinstance(stretch, LinearStretch)
            and stretch.slope == 1 and stretch.intercept == 0):
        return ColorDist.LinearDist(hashsize)
    return _AstropyStretchDist(hashsize, stretch)


def docs_from_super_if_missing(cls):
    """
    Decorator to copy the docstrings from the interface methods to the
    methods in the class.
    """
    for name, method in cls.__dict__.items():
        if not name.startswith("_"):
            if method.__doc__:
                continue
            interface_method = getattr(ImageViewerLogic, name, None)
            if interface_method:
                method.__doc__ = interface_method.__doc__
    return cls


# The inheritance order below matters -- VBox needs to come first
@docs_from_super_if_missing
class ImageWidget(ipyw.VBox, ImageViewerLogic):
    """
    Image widget for Jupyter notebook using the Ginga viewer.

    This widget implements the astro-image-display-api (AIDA)
    `~astro_image_display_api.ImageViewerInterface`, plus a handful of
    interactive conveniences that Ginga supports natively (live cursor
    readout, click-to-center, and interactive marking).

    Parameters
    ----------
    logger : obj or ``None``
        Ginga logger. For example::

            from ginga.misc.log import get_logger
            logger = get_logger('my_viewer', log_stderr=False,
                                log_file='ginga.log', level=40)

    image_width, image_height : int
        Dimension of the Jupyter notebook's image widget, in pixels.
    """
    # Allowed locations for cursor display
    ALLOWED_CURSOR_LOCATIONS = ['top', 'bottom', None]

    # List of marker names that are for internal use only
    RESERVED_MARKER_SET_NAMES = ['all']

    # Default marker name for marking via API
    DEFAULT_MARKER_NAME: str = 'default'

    # Default marker name for interactive marking
    DEFAULT_INTERACTIVE_MARKER_NAME: str = 'interactive'

    def __init__(self, *args, logger=None, image_width=500, image_height=500,
                 **kwargs):
        if 'use_opencv' in kwargs:
            # Pop it so it is not forwarded to VBox, which would raise a
            # TraitError for an unknown trait.
            kwargs.pop('use_opencv')
            warnings.warn("use_opencv kwarg has been deprecated--"
                          "opencv will be used if it is installed",
                          DeprecationWarning)

        super().__init__(*args, **kwargs)
        # ImageViewerLogic is a dataclass; we do not run its __init__, so run
        # the post-init hook it would otherwise provide to set up its state.
        ImageViewerLogic.__post_init__(self)

        # Ginga displays a single image at a time; remember which label is on
        # screen so per-label settings for other images do not touch the
        # display (see _is_displayed_label).
        self._displayed_image_label = None

        self._viewer = EnhancedCanvasView(logger=logger)

        self._jup_img = ipyw.Image(format='jpeg')

        # Set the image margin over the widget's default of 2px on all sides.
        self._jup_img.layout.margin = '0'

        # Set both of these to ensure consistent display in the notebook and in
        # jupyterlab when the image is put into a container smaller than itself.
        self._jup_img.max_width = '100%'
        self._jup_img.height = 'auto'

        # Set the width of the box containing the image to the desired width.
        # The height is set automatically by the image aspect ratio.
        self.layout.width = str(image_width)

        # Ginga uses these to figure out what size image to make.
        self._jup_img.width = image_width
        self._jup_img.height = image_height

        self._viewer.set_widget(self._jup_img)

        # enable all possible keyboard and pointer operations
        self._viewer.get_bindings().enable_all(True)

        # Shapes used to draw catalog markers on the viewer's canvas.
        self.dc = drawCatalog

        # Internal interaction-state trackers; all start disabled.
        self._is_marking = False
        self._click_center = False
        self._click_drag = False
        self._scroll_pan = False

        # Match the ginga defaults.
        self.scroll_pan = True
        self.click_drag = False

        bind_map = self._viewer.get_bindmap()
        # Right-click and drag adjusts the contrast; shift-right-click restores.
        bind_map.map_event(None, (), 'ms_right', 'contrast')
        bind_map.map_event(None, ('shift',), 'ms_right', 'contrast_restore')

        # State for interactive marking.
        self._interactive_marker_set_name = self.DEFAULT_INTERACTIVE_MARKER_NAME
        self._interactive_points = []
        self._interactive_style = self._default_catalog_style.copy()

        # Guards re-entrancy while we programmatically update the viewport so
        # that the live-state sync in get_viewport does not clobber the values
        # we just set.
        self._updating_viewport = False

        # coordinates display
        self._jup_coord = ipyw.HTML('Coordinates show up here')
        self._viewer.add_callback('cursor-changed', self._mouse_move_cb)
        self._viewer.add_callback('cursor-down', self._mouse_click_cb)

        # Output widget that captures printed output, for debugging.
        self._print_out = ipyw.Output()

        self._cursor = 'bottom'
        self.children = [self._jup_img, self._jup_coord]

    # ------------------------------------------------------------------
    # Read-only conveniences
    # ------------------------------------------------------------------
    @property
    def logger(self):
        """Logger for this widget."""
        return self._viewer.logger

    @property
    def viewer(self):
        """The underlying ginga viewer object."""
        return self._viewer

    @property
    def image_width(self):
        """Width of the image widget, in pixels."""
        return int(self._jup_img.width)

    @property
    def image_height(self):
        """Height of the image widget, in pixels."""
        return int(self._jup_img.height)

    @property
    def print_out(self):
        """
        An `ipywidgets.Output` widget that captures printed output from the
        viewer, intended primarily for debugging.
        """
        return self._print_out

    # ------------------------------------------------------------------
    # Mouse callbacks
    # ------------------------------------------------------------------
    def _mouse_move_cb(self, viewer, button, data_x, data_y):
        """
        Callback to display the cursor position (and value) in the readout.
        """
        if self.cursor is None:  # no-op
            return

        image = viewer.get_image()
        if image is not None:
            ix = int(data_x + 0.5)
            iy = int(data_y + 0.5)
            try:
                imval = viewer.get_data(ix, iy)
                imval = '{:8.3f}'.format(imval)
            except Exception:
                imval = 'N/A'

            val = 'X: {:.2f}, Y: {:.2f}'.format(data_x, data_y)
            if image.wcs.wcs is not None:
                try:
                    ra, dec = image.pixtoradec(data_x, data_y)
                    val += ' (RA: {}, DEC: {})'.format(
                        ra_deg_to_str(ra), dec_deg_to_str(dec))
                except Exception:
                    val += ' (RA, DEC: WCS error)'

            val += ', value: {}'.format(imval)
            self._jup_coord.value = val

    def _mouse_click_cb(self, viewer, event, data_x, data_y):
        """
        Callback to handle mouse clicks for marking or centering.
        """
        if self.is_marking:
            self._append_interactive_marker(data_x, data_y)
            with self._print_out:
                print('Selected {} {}'.format(data_x, data_y))
        elif self.click_center:
            self._center_on((data_x, data_y))
            with self._print_out:
                print('Centered on X={} Y={}'.format(data_x, data_y))

    # ------------------------------------------------------------------
    # Image loading
    # ------------------------------------------------------------------
    def load_image(self, image, image_label=None, **kwargs):
        # Let the AIDA logic store the data + WCS and set up the initial
        # viewport, cuts and stretch in our internal state.
        super().load_image(image, image_label=image_label, **kwargs)

        # The just-loaded image is now the one on screen.
        self._displayed_image_label = self._resolve_image_label(image_label)

        # Build a ginga AstroImage from the stored data and show it. The
        # helpers resolve image_label internally, matching the bqplot backend.
        data = self.get_image(image_label=image_label)
        self._viewer.set_image(self._build_ginga_image(data))

        # Apply the stored viewport, cuts, stretch and colormap to the freshly
        # displayed image.
        self._apply_viewport_to_ginga(image_label)
        self._apply_cuts_to_ginga(image_label)
        self._apply_stretch_to_ginga(image_label)
        self._apply_colormap_to_ginga(image_label)

    def _is_displayed_label(self, image_label):
        """
        Return `True` if ``image_label`` refers to the image currently shown
        in the ginga viewer.

        The ginga viewer shows one image at a time, but the AIDA state stores
        settings per label. Helpers that push stored state into the live
        viewer must be skipped when the label refers to a different,
        not-displayed image; its stored settings are applied when it is next
        loaded.
        """
        try:
            image_label = self._resolve_image_label(image_label)
        except ValueError:
            return False
        return image_label == self._displayed_image_label

    def _build_ginga_image(self, data):
        """
        Build a ginga `~ginga.AstroImage.AstroImage` from stored image data,
        carrying over the WCS when one is available.
        """
        image = AstroImage(logger=self.logger)

        if isinstance(data, NDData):
            image.set_data(np.asarray(data.data))
            wcs = getattr(data, 'wcs', None)
            if wcs is not None:
                # A bad or exotic WCS can fail anywhere in here -- header
                # serialization, ginga's header parser, or set_wcs -- so guard
                # the whole block and degrade gracefully to displaying the image
                # without a WCS. Catch Exception (not BaseException) so Ctrl-C
                # still propagates, and log the traceback for debugging.
                try:
                    from ginga.util.wcsmod.wcs_astropy import AstropyWCS
                    _wcs = AstropyWCS(self.logger)
                    _wcs.load_header(wcs.to_header())
                    image.set_wcs(_wcs)
                except Exception:  # pragma: no cover - defensive
                    self.logger.warning('Unable to set WCS from image',
                                        exc_info=True)
        else:
            image.set_data(np.asarray(data))

        return image

    # ------------------------------------------------------------------
    # Cuts / stretch / colormap
    # ------------------------------------------------------------------
    def set_cuts(self, value, image_label=None, **kwargs):
        super().set_cuts(value, image_label=image_label, **kwargs)
        # Changing the cuts only affects the color mapping, so leave the
        # current viewport (zoom/pan) untouched.
        self._apply_cuts_to_ginga(image_label)

    def _apply_cuts_to_ginga(self, image_label):
        ginga_image = self._viewer.get_image()
        if ginga_image is None or not self._is_displayed_label(image_label):
            return
        cuts = self.get_cuts(image_label=image_label)
        # get_cuts always returns a BaseInterval, and get_limits computes the
        # concrete low/high for any interval type (percentile, zscale, etc.),
        # so this is exact even for non-min/max cuts.
        low, high = cuts.get_limits(ginga_image.get_data())
        self._viewer.cut_levels(low, high)

    def set_stretch(self, value, image_label=None, **kwargs):
        # The astropy stretch (including its shape parameter) is reproduced in
        # ginga by configuring the matching native color distribution, or, for
        # stretches ginga lacks a family for, by an astropy-backed adapter.
        # See _ginga_dist_for_stretch.
        super().set_stretch(value, image_label=image_label, **kwargs)
        # Changing the stretch only affects the color mapping, so leave the
        # current viewport (zoom/pan) untouched.
        self._apply_stretch_to_ginga(image_label)

    def _apply_stretch_to_ginga(self, image_label):
        if (self._viewer.get_image() is None
                or not self._is_displayed_label(image_label)):
            return
        stretch = self.get_stretch(image_label=image_label)
        rgbmap = self._viewer.get_rgbmap()
        # Inject a fully parametrized distribution; set_dist fires the rgbmap's
        # 'changed' callback, which redraws. (set_color_algorithm cannot carry
        # parameters -- it only ever builds a distribution with ginga defaults.)
        rgbmap.set_dist(_ginga_dist_for_stretch(stretch, rgbmap.get_hash_size()))

    def set_colormap(self, map_name, image_label=None, **kwargs):
        # Ginga only logs (rather than raises) on an unknown colormap name,
        # which would leave the stored state disagreeing with the display, so
        # validate the name up front.
        if map_name not in ginga_cmap.get_names():
            raise ValueError(f'Colormap {map_name!r} is not a valid ginga '
                             'colormap name.')
        super().set_colormap(map_name, image_label=image_label, **kwargs)
        self._apply_colormap_to_ginga(image_label)

    def _apply_colormap_to_ginga(self, image_label):
        if (self._viewer.get_image() is None
                or not self._is_displayed_label(image_label)):
            return
        cmap_name = self.get_colormap(image_label=image_label)
        # The colormap defaults to None and is not set on load, so guard
        # against pushing a None into ginga's set_color_map.
        if cmap_name is not None:
            self._viewer.set_color_map(cmap_name)

    # ------------------------------------------------------------------
    # Catalog API
    # ------------------------------------------------------------------
    def load_catalog(self, table, **kwargs):
        super().load_catalog(table, **kwargs)
        catalog_label = kwargs.get("catalog_label", None)
        self._draw_catalog(catalog_label)

    def set_catalog_style(self, catalog_label=None, shape="circle",
                          color="red", size=5, **kwargs):
        super().set_catalog_style(
            catalog_label=catalog_label, shape=shape, color=color, size=size,
            **kwargs
        )
        self._draw_catalog(catalog_label)

    def remove_catalog(self, catalog_label=None, **kwargs):
        # Figure out which canvas tags to remove before the bookkeeping in
        # super() forgets about them.
        if catalog_label == "*":
            tags = [str(label) for label in self._catalogs]
        else:
            tags = [str(self._resolve_catalog_label(catalog_label))]

        super().remove_catalog(catalog_label, **kwargs)

        for tag in tags:
            try:
                self._viewer.canvas.delete_object_by_tag(tag)
            except Exception:
                # Nothing was drawn under this tag; that is fine.
                pass

    def _make_marker(self, style):
        """
        Build a ginga drawing-object factory from a catalog style dict.
        """
        shape = style.get('shape', 'circle')
        color = style.get('color', 'red')
        size = style.get('size', 5)
        linewidth = style.get('linewidth', 1)

        if shape in ('square', 'box'):
            return functools.partial(self.dc.SquareBox, radius=size, color=color,
                                     linewidth=linewidth)
        elif shape in ('cross', 'plus', 'crosshair'):
            point_style = 'plus' if shape == 'plus' else 'cross'
            return functools.partial(self.dc.Point, radius=size, style=point_style,
                                     color=color, linewidth=linewidth)
        else:  # circle and anything we do not specifically handle
            return functools.partial(self.dc.Circle, radius=size, color=color,
                                     linewidth=linewidth)

    def _draw_catalog(self, catalog_label):
        """
        Draw (or redraw) a catalog's markers on the ginga canvas under a tag
        derived from the catalog label.
        """
        tag = str(self._resolve_catalog_label(catalog_label))

        # Remove any existing drawing for this catalog before redrawing.
        try:
            self._viewer.canvas.delete_object_by_tag(tag)
        except Exception:
            pass

        catalog = self.get_catalog(catalog_label=catalog_label)
        if catalog is None or len(catalog) == 0:
            return

        marker = self._make_marker(self.get_catalog_style(catalog_label=catalog_label))

        objs = []
        for x, y in zip(catalog['x'], catalog['y']):
            if x is None or y is None or np.ma.is_masked(x) or np.ma.is_masked(y):
                continue
            objs.append(marker(x=float(x), y=float(y), coord='data'))

        if objs:
            self._viewer.canvas.add(self.dc.CompoundObject(*objs), tag=tag)

    # ------------------------------------------------------------------
    # Viewport API
    # ------------------------------------------------------------------
    @property
    def _viewport_window_size(self):
        """
        The window dimension (in screen pixels) used to convert between a
        pixel field of view and a ginga scale.

        Ginga scales the image isotropically (the same screen-pixels-per-data-
        pixel in both directions), so we define the field of view as the
        *horizontal* extent and always use the window width as the reference
        dimension. This is self-consistent for a non-square viewer: both
        set_viewport and get_viewport go through this single dimension, so the
        round-trip is exact regardless of the window aspect ratio.
        """
        return self._viewer.get_window_size()[0]

    def set_viewport(self, center=None, fov=None, image_label=None, **kwargs):
        # The AIDA logic handles all of the WCS bookkeeping (and validation),
        # which we need in case center/fov are in sky coordinates.
        super().set_viewport(center=center, fov=fov, image_label=image_label,
                             **kwargs)
        self._apply_viewport_to_ginga(image_label)

    def _apply_viewport_to_ginga(self, image_label):
        """
        Push the stored viewport (center + fov) into the ginga viewer as a pan
        position and scale.
        """
        if (self._viewer.get_image() is None
                or not self._is_displayed_label(image_label)):
            return

        # Read the stored viewport back in pixel coordinates without going
        # through our own get_viewport (which would try to sync ginga state).
        viewport = super().get_viewport(sky_or_pixel='pixel',
                                        image_label=image_label)
        center_x, center_y = viewport['center']
        fov_pixels = viewport['fov']

        self._updating_viewport = True
        try:
            self._viewer.set_pan(center_x, center_y)
            if fov_pixels:
                scale = self._viewport_window_size / fov_pixels
                self._viewer.scale_to(scale, scale)
        finally:
            self._updating_viewport = False

    def get_viewport(self, sky_or_pixel=None, image_label=None, **kwargs):
        # Reflect any interactive pan/zoom the user did in the browser by
        # syncing the live ginga state into our stored viewport before
        # returning it. This mirrors the bqplot backend, which reflects
        # interactive viewport changes in get_viewport.
        if not self._updating_viewport:
            self._sync_ginga_to_stored_viewport(image_label)
        return super().get_viewport(sky_or_pixel, image_label, **kwargs)

    def _sync_ginga_to_stored_viewport(self, image_label):
        """
        Update the stored viewport from the live ginga pan/scale, but only if
        the user has actually panned or zoomed since we last set it. Skipping
        the update when nothing has changed keeps programmatically-set values
        exact (important for sky/pixel round-trips).
        """
        try:
            image_label = self._resolve_image_label(image_label)
        except ValueError:
            # Ambiguous or unknown label -- let super().get_viewport raise the
            # appropriate error instead of masking it here.
            return

        if (image_label not in self._images
                or self._viewer.get_image() is None
                or image_label != self._displayed_image_label):
            # Only the displayed image can have been interactively panned or
            # zoomed; the live ginga state belongs to it alone.
            return

        # What pan/scale does the currently-stored viewport correspond to?
        stored = super().get_viewport(sky_or_pixel='pixel',
                                      image_label=image_label)
        exp_x, exp_y = stored['center']
        exp_fov = stored['fov']
        window = self._viewport_window_size
        exp_scale = window / exp_fov if exp_fov else None

        cur_x, cur_y = self._viewer.get_pan()
        cur_scale = self._viewer.get_scale()

        changed = (
            abs(cur_x - exp_x) > 1e-6
            or abs(cur_y - exp_y) > 1e-6
            or (exp_scale is not None and abs(cur_scale - exp_scale) > 1e-9)
        )
        if not changed:
            return

        # The user interacted; store the new pixel center and field of view.
        fov_pixels = window / cur_scale if cur_scale else exp_fov
        self._updating_viewport = True
        try:
            super().set_viewport(center=(cur_x, cur_y), fov=fov_pixels,
                                 image_label=image_label)
        finally:
            self._updating_viewport = False

    # ------------------------------------------------------------------
    # Interactive extras: cursor, click-to-center, marking
    # ------------------------------------------------------------------
    def _center_on(self, point):
        """
        Center the view on a point given either as a `~astropy.coordinates.SkyCoord`
        or as a tuple of pixel ``(X, Y)`` coordinates.
        """
        self.set_viewport(center=point)

    @property
    def cursor(self):
        """
        Show or hide cursor information (X, Y, WCS).
        Acceptable values are 'top', 'bottom', or ``None``.
        """
        return self._cursor

    @cursor.setter
    def cursor(self, val):
        if val is None:
            self._jup_coord.layout.visibility = 'hidden'
            self._jup_coord.layout.display = 'none'
        elif val == 'top' or val == 'bottom':
            self._jup_coord.layout.visibility = 'visible'
            self._jup_coord.layout.display = 'flex'
            if val == 'top':
                self.layout.flex_flow = 'column-reverse'
            else:
                self.layout.flex_flow = 'column'
        else:
            raise ValueError('Invalid value {} for cursor. Valid values are: '
                             '{}'.format(val, self.ALLOWED_CURSOR_LOCATIONS))
        self._cursor = val

    @property
    def click_center(self):
        """
        Settable. If `True`, clicking on the image centers the view on the
        clicked point. If `False`, that interaction is disabled.
        """
        return self._click_center

    @click_center.setter
    def click_center(self, val):
        if not isinstance(val, bool):
            raise ValueError('Must be True or False')
        elif self.is_marking and val:
            raise ValueError('Cannot set to True while marking is active')

        if val:
            self.click_drag = False

        self._click_center = val

    @property
    def click_drag(self):
        """
        Settable. If `True`, click-and-drag pans the image. If `False`, it is
        disabled. Automatically made `False` while marking is active.
        """
        return self._click_drag

    @click_drag.setter
    def click_drag(self, value):
        if not isinstance(value, bool):
            raise ValueError('click_drag must be either True or False')
        if self.is_marking:
            raise ValueError('Interactive marking is in progress. Call '
                             'stop_marking() to end marking before setting '
                             'click_drag')
        self._click_drag = value
        bindmap = self._viewer.get_bindmap()
        if value:
            # Only turn off click_center if click_drag is being set to True.
            self.click_center = False
            bindmap.map_event(None, (), 'ms_left', 'pan')
        else:
            bindmap.map_event(None, (), 'ms_left', 'cursor')

    @property
    def scroll_pan(self):
        """
        Settable. If `True`, scrolling pans the image. If `False`, scrolling
        (up/down) zooms the image in and out.
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

    @property
    def is_marking(self):
        """
        `True` if in interactive marking mode, `False` otherwise. While
        marking, a mouse click adds a new marker to the interactive catalog.
        """
        return self._is_marking

    def start_marking(self, marker_name=None, marker=None):
        """
        Begin interactive marking.

        Parameters
        ----------
        marker_name : str, optional
            Catalog label under which the interactively-placed markers are
            stored. Retrieve them later with ``get_catalog(catalog_label=...)``.

        marker : dict, optional
            Style for the interactive markers, e.g.
            ``{'shape': 'circle', 'color': 'cyan', 'size': 20}``.
        """
        self._cached_state = dict(click_center=self.click_center,
                                  click_drag=self.click_drag,
                                  scroll_pan=self.scroll_pan)
        self.click_center = False
        self.click_drag = False
        # Ensure there is still a mouse way to pan.
        self.scroll_pan = True
        self._is_marking = True
        self._interactive_points = []

        if marker_name is not None:
            self._validate_marker_name(marker_name)
            self._interactive_marker_set_name = marker_name
        else:
            self._interactive_marker_set_name = self.DEFAULT_INTERACTIVE_MARKER_NAME

        if marker is not None:
            self._interactive_style = marker

    def stop_marking(self, clear_markers=False):
        """
        Stop interactive marking.

        Parameters
        ----------
        clear_markers : bool, optional
            If `True`, remove the interactively-placed markers. If `False`
            (default), they are retained and can be retrieved with
            ``get_catalog``.
        """
        if self.is_marking:
            self._is_marking = False
            self.click_center = self._cached_state['click_center']
            self.click_drag = self._cached_state['click_drag']
            self.scroll_pan = self._cached_state['scroll_pan']
            self._cached_state = {}
            if clear_markers:
                try:
                    self.remove_catalog(
                        catalog_label=self._interactive_marker_set_name)
                except ValueError:
                    pass
                self._interactive_points = []

    def _validate_marker_name(self, marker_name):
        """Raise an error if ``marker_name`` is reserved."""
        if marker_name in self.RESERVED_MARKER_SET_NAMES:
            raise ValueError('The marker name {} is not allowed. Any name is '
                             'allowed except these: '
                             '{}'.format(marker_name,
                                         ', '.join(self.RESERVED_MARKER_SET_NAMES)))

    def _append_interactive_marker(self, x, y):
        """Add a point to the interactive catalog and redraw it."""
        self._interactive_points.append((x, y))
        table = Table(rows=self._interactive_points, names=['x', 'y'])
        self.load_catalog(table,
                          catalog_label=self._interactive_marker_set_name,
                          catalog_style=self._interactive_style)

    # ------------------------------------------------------------------
    # Saving
    # ------------------------------------------------------------------
    def save(self, filename, overwrite=False, **kwargs):
        """
        Save the current image view to the given image file.

        Parameters
        ----------
        filename : str or `os.PathLike`
            Name of the file to save to.

        overwrite : bool, optional
            If `True`, overwrite an existing file.
        """
        if not overwrite and Path(filename).exists():
            raise FileExistsError(f'File {filename} exists and overwrite=False')

        # Ginga renders the view server-side, so this works without a running
        # browser and honors the image format implied by the file extension.
        self._viewer.save_rgb_image_as_file(str(filename))
