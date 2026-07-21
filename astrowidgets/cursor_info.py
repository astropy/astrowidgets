"""
Shared cursor-information readout for astrowidgets backends.

Backends compose their image widget with an HTML readout that shows the
cursor position (X/Y and, when a WCS is available, RA/Dec) and the image
value under the cursor. The formatting and widget plumbing live here so
that every backend behaves the same; each backend only forwards its
native mouse-move events to `CursorInfoMixin._update_cursor_text`.
"""
import numpy as np

from astropy import units as u
from astropy.nddata import NDData

import ipywidgets as ipw

__all__ = ['format_cursor_text', 'CursorInfoMixin', 'READOUT_TEMPLATE',
           'ALLOWED_CURSOR_LOCATIONS', 'ALLOWED_SKY_COORDINATE_FORMATS']

ALLOWED_CURSOR_LOCATIONS = ('top', 'bottom', None)
ALLOWED_SKY_COORDINATE_FORMATS = ('degrees', 'sexagesimal')

# The readout is rendered in a <pre> so that the fixed-width padding in
# format_cursor_text survives HTML whitespace collapsing and the digits
# line up, keeping the line from jittering as the cursor moves.
READOUT_TEMPLATE = '<pre style="margin: 0">{}</pre>'


def format_cursor_text(x, y, data=None, wcs=None, sky_format='degrees'):
    """
    Format the cursor readout for a position in an image.

    Parameters
    ----------
    x, y : float
        Cursor position in data (pixel) coordinates, with pixel centers
        at integer values (origin 0), i.e. the convention of
        `astropy.wcs.WCS.pixel_to_world`.
    data : array-like or None
        Image data, indexed as ``data[y, x]``. When ``None``, or when the
        cursor is outside the array, the value reads ``N/A``.
    wcs : `astropy.wcs.WCS` or None
        WCS used for the RA/Dec segment; omitted when ``None``. The sky
        position is always converted to ICRS, whatever the native frame
        of the WCS, and the segment is tagged ``(ICRS)`` to say so.
    sky_format : str
        One of `ALLOWED_SKY_COORDINATE_FORMATS`. ``'degrees'`` shows
        decimal degrees, ``'sexagesimal'`` shows ``HH:MM:SS.ss`` /
        ``±DD:MM:SS.ss``.

    Returns
    -------
    str
        Readout like ``'X: 10    Y: 7     RA: 137.0412 Dec:
        -10.1235 (ICRS) value: 3.1'``. X/Y show the integer pixel the
        value is sampled from. The X/Y and decimal RA/Dec fields are
        fixed width, left-aligned so the padding trails the numbers and
        the line does not jitter as the cursor moves (sexagesimal is
        naturally fixed width); at their widest, fields are followed by
        a single space.
    """
    if sky_format not in ALLOWED_SKY_COORDINATE_FORMATS:
        raise ValueError(f'Invalid value {sky_format!r} for sky_format. '
                         f'Valid values are: {ALLOWED_SKY_COORDINATE_FORMATS}')

    if data is not None:
        data = np.asarray(data)

    # floor rather than int() so that positions just outside the image on
    # the negative side do not alias onto pixel 0.
    x_index = int(np.floor(x + 0.5))
    y_index = int(np.floor(y + 0.5))

    if (data is not None
            and 0 <= y_index < data.shape[0]
            and 0 <= x_index < data.shape[1]):
        value = f'value: {data[y_index, x_index]:.1f}'
    else:
        value = 'value: N/A'

    segments = [f'X: {x_index:<5d} Y: {y_index:<5d}']

    if wcs is not None:
        try:
            sky = wcs.pixel_to_world(x, y).icrs
            if sky_format == 'sexagesimal':
                ra = sky.ra.to_string(unit=u.hour, sep=':',
                                      precision=2, pad=True)
                dec = sky.dec.to_string(unit=u.degree, sep=':',
                                        precision=2, alwayssign=True,
                                        pad=True)
            else:
                # Fixed width (RA 0-360, Dec always signed) so the line
                # does not jitter as the cursor moves.
                ra = f'{sky.ra.deg:<8.4f}'
                dec = f'{sky.dec.deg:<+8.4f}'
            segments.append(f'RA: {ra} Dec: {dec} (ICRS)')
        except Exception:
            # Deliberately broad: the WCS is user-supplied and only
            # duck-typed to the APE-14 interface, so the exceptions it can
            # raise are open-ended, and this runs on every mouse move,
            # where an escaped exception would spam the log or vanish
            # instead of surfacing cleanly. Degrade to a visible error.
            segments.append('RA/Dec: WCS error')

    segments.append(value)
    return ' '.join(segments)


class CursorInfoMixin:
    """
    Cursor-information readout shared by the astrowidgets backends.

    Intended for `ipywidgets.VBox` subclasses that also inherit
    `~astro_image_display_api.image_viewer_logic.ImageViewerLogic`, whose
    per-label image storage supplies the data and WCS for the readout.
    The mixin defines no ``__init__``; backends call `_init_cursor_info`
    from their own ``__init__``, place the returned readout widget in
    their own ``children``, and forward their native mouse-move events
    to `_update_cursor_text`.
    """
    ALLOWED_CURSOR_LOCATIONS = ALLOWED_CURSOR_LOCATIONS
    ALLOWED_SKY_COORDINATE_FORMATS = ALLOWED_SKY_COORDINATE_FORMATS

    def _init_cursor_info(self):
        """
        Create and return the readout widget. The backend adds it to its
        own ``children``, below the image widget; the `cursor` property
        flips the box to ``column-reverse`` to show it on top.
        """
        self._cursor_readout = ipw.HTML('Coordinates show up here')
        self._sky_coordinate_format = 'degrees'
        self._last_cursor_position = None
        self.cursor = 'bottom'
        return self._cursor_readout

    @property
    def cursor(self):
        """
        Show or hide cursor information (X, Y, RA/Dec, value).
        Acceptable values are 'top', 'bottom', or ``None``.
        """
        return self._cursor_location

    @cursor.setter
    def cursor(self, val):
        """
        Set the readout location.

        Parameters
        ----------
        val : str or None
            One of `ALLOWED_CURSOR_LOCATIONS`: ``'top'`` or ``'bottom'``
            place the readout above or below the image by flipping the
            box's flex direction; ``None`` hides it. The readout is
            re-rendered at the last cursor position. Anything else
            raises `ValueError`.
        """
        if val is None:
            self._cursor_readout.layout.visibility = 'hidden'
            self._cursor_readout.layout.display = 'none'
        elif val == 'top' or val == 'bottom':
            self._cursor_readout.layout.visibility = 'visible'
            self._cursor_readout.layout.display = 'flex'
            if val == 'top':
                self.layout.flex_flow = 'column-reverse'
            else:
                self.layout.flex_flow = 'column'
        else:
            raise ValueError('Invalid value {} for cursor. Valid values are: '
                             '{}'.format(val, self.ALLOWED_CURSOR_LOCATIONS))
        self._cursor_location = val
        self._refresh_cursor_text()

    @property
    def sky_coordinate_format(self):
        """
        Format of the RA/Dec part of the cursor readout, either
        ``'degrees'`` (decimal degrees, the default) or ``'sexagesimal'``.
        """
        return self._sky_coordinate_format

    @sky_coordinate_format.setter
    def sky_coordinate_format(self, val):
        """
        Set the RA/Dec display format.

        Parameters
        ----------
        val : str
            One of `ALLOWED_SKY_COORDINATE_FORMATS`: ``'degrees'`` or
            ``'sexagesimal'``. The readout is immediately re-rendered
            at the last cursor position, so the change is visible
            without waiting for the next mouse move. Anything else
            raises `ValueError`.
        """
        if val not in self.ALLOWED_SKY_COORDINATE_FORMATS:
            raise ValueError('Invalid value {} for sky_coordinate_format. '
                             'Valid values are: '
                             '{}'.format(val,
                                         self.ALLOWED_SKY_COORDINATE_FORMATS))
        self._sky_coordinate_format = val
        self._refresh_cursor_text()

    def _update_cursor_text(self, x, y):
        """
        Update the readout for a cursor at ``(x, y)`` in data coordinates,
        using the data and WCS of the currently displayed image.

        Parameters
        ----------
        x, y : float
            Cursor position in data (pixel) coordinates, following the
            `astropy.wcs.WCS.pixel_to_world` convention (pixel centers
            at integer values, origin 0).
        """
        self._last_cursor_position = (x, y)
        self._refresh_cursor_text()

    def _refresh_cursor_text(self):
        """
        Re-render the readout at the last known cursor position, e.g.
        after the sky coordinate format changes.
        """
        if (self._last_cursor_position is None or self.cursor is None
                or not self._displayed_image_labels):
            return

        x, y = self._last_cursor_position
        info = self._images[self._displayed_image_labels[0]]
        data = info.data
        if isinstance(data, NDData):
            data = data.data

        self._cursor_readout.value = READOUT_TEMPLATE.format(
            format_cursor_text(x, y, data=data, wcs=info.wcs,
                               sky_format=self._sky_coordinate_format))
