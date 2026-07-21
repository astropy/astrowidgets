import numpy as np
import pytest

from astropy import units as u
from astropy.nddata import NDData
from astropy.wcs import WCS

from astrowidgets.cursor_info import format_cursor_text, READOUT_TEMPLATE


# The wcs fixture these tests use is shared across test modules; it lives
# in conftest.py.

@pytest.fixture
def data():
    # Non-square so that a transposed x/y lookup is detectable.
    return np.arange(35, dtype=float).reshape(5, 7)


def test_format_pixel_and_value(data):
    text = format_cursor_text(3.0, 2.0, data=data)
    # data[y, x] = data[2, 3] = 17. Exact string: X/Y are width-5
    # left-aligned fields with a single space of separation, so a
    # 5-digit pixel index is followed by exactly one space.
    assert text == 'X: 3     Y: 2     value: 17.0'


def test_format_uses_correct_pixel_to_world_order(wcs, data):
    # Regression test: the bqplot backend used to call
    # wcs.pixel_to_world(y, x), transposing the sky coordinates.
    x, y = 3.0, 2.0
    good = wcs.pixel_to_world(x, y).icrs
    swapped = wcs.pixel_to_world(y, x).icrs
    assert good.ra.deg != swapped.ra.deg

    text = format_cursor_text(x, y, data=data, wcs=wcs)
    assert f'RA: {good.ra.deg:<8.4f}' in text
    assert f'Dec: {good.dec.deg:<+8.4f}' in text


def test_format_sexagesimal(wcs):
    # The sexagesimal readout must match astropy's own to_string
    # rendering: RA as HH:MM:SS.ss (hours, zero-padded), Dec as
    # +/-DD:MM:SS.ss (always signed).
    x, y = 3.0, 2.0
    sky = wcs.pixel_to_world(x, y).icrs
    text = format_cursor_text(x, y, wcs=wcs, sky_format='sexagesimal')
    ra = sky.ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
    dec = sky.dec.to_string(unit=u.degree, sep=':', precision=2,
                            alwayssign=True, pad=True)
    assert f'RA: {ra}' in text
    assert f'Dec: {dec}' in text


def test_format_galactic_wcs_converts_to_icrs_and_says_so():
    # The readout always shows ICRS, even for a galactic-native WCS
    # (e.g. the Spitzer example image), and tags the segment to say so.
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [50, 50]
    wcs.wcs.cdelt = np.array([-0.001, 0.001])
    wcs.wcs.crval = [18.35, 0.16]
    wcs.wcs.ctype = ['GLON-CAR', 'GLAT-CAR']

    x, y = 3.0, 2.0
    sky = wcs.pixel_to_world(x, y)
    assert sky.frame.name == 'galactic'
    icrs = sky.icrs

    text = format_cursor_text(x, y, wcs=wcs)
    assert f'RA: {icrs.ra.deg:<8.4f}' in text
    assert f'Dec: {icrs.dec.deg:<+8.4f} (ICRS)' in text


def test_format_invalid_sky_format():
    # An unrecognized sky_format must raise rather than silently fall
    # back to a default, and the error must name the bad argument.
    with pytest.raises(ValueError, match='sky_format'):
        format_cursor_text(1.0, 1.0, sky_format='hours')


@pytest.mark.parametrize('x,y', [(7.0, 2.0), (3.0, 5.0), (-0.6, 2.0)])
def test_format_out_of_bounds_value(x, y, data):
    text = format_cursor_text(x, y, data=data)
    assert 'value: N/A' in text


def test_format_decimal_readout_is_fixed_width(wcs):
    # The X/Y and decimal RA/Dec fields are padded to a constant width
    # so the readout line does not jitter side-to-side as the cursor
    # moves, including across changes in digit count (9.99 -> 10.01).
    big_data = np.zeros((200, 300))
    lengths = {len(format_cursor_text(x, 2.0, data=big_data, wcs=wcs))
               for x in (0.0, 9.99, 10.01, 100.0)}
    assert len(lengths) == 1


def test_format_displayed_pixel_matches_sampled_pixel(data):
    # The X/Y readout shows the integer pixel the value is sampled from,
    # not the fractional cursor position.
    text = format_cursor_text(2.6, 1.4, data=data)
    assert f'X: {3:<5d}' in text
    assert f'Y: {1:<5d}' in text
    # data[y, x] = data[1, 3] = 10
    assert 'value: 10.0' in text


def test_format_negative_coordinate_does_not_alias_to_pixel_zero(data):
    # int(x + 0.5) truncates toward zero, so a cursor at x=-0.4 would
    # wrongly read pixel 0; floor-based rounding must report N/A instead.
    text = format_cursor_text(-0.6, 0.0, data=data)
    assert 'value: N/A' in text


def test_format_accepts_list_data():
    # data is documented as array-like, so a plain nested list must work,
    # not just objects that already have .shape and 2D indexing.
    text = format_cursor_text(1.0, 1.0, data=[[0.0, 1.0], [2.0, 3.0]])
    assert 'value: 3.0' in text


def test_format_no_data():
    text = format_cursor_text(1.0, 1.0)
    assert 'value: N/A' in text


def test_format_wcs_error(data):
    class BrokenWCS:
        def pixel_to_world(self, x, y):
            raise ValueError('broken')

    text = format_cursor_text(1.0, 1.0, data=data, wcs=BrokenWCS())
    assert 'RA/Dec: WCS error' in text


# ----------------------------------------------------------------------
# Widget behavior, shared across backends
# ----------------------------------------------------------------------

def _backend_widget(backend):
    if backend == 'bqplot':
        pytest.importorskip('bqplot')
        from astrowidgets.bqplot import ImageWidget
    else:
        pytest.importorskip('ginga')
        from astrowidgets.ginga import ImageWidget
    return ImageWidget()


def _simulate_mouse_move(widget, backend, x, y):
    if backend == 'bqplot':
        # Simulate the comm message a real mouse move produces.
        widget._astro_im.interaction._handle_custom_msg(
            {'event': 'mousemove', 'domain': {'x': x, 'y': y}}, [])
    else:
        widget._mouse_move_cb(widget._viewer, None, x, y)


@pytest.mark.parametrize('backend', ['bqplot', 'ginga'])
def test_cursor_property(backend):
    widget = _backend_widget(backend)
    assert widget.cursor == 'bottom'
    assert widget.layout.flex_flow == 'column'

    widget.cursor = 'top'
    assert widget.layout.flex_flow == 'column-reverse'
    assert widget._cursor_readout.layout.visibility == 'visible'

    widget.cursor = None
    assert widget._cursor_readout.layout.display == 'none'

    with pytest.raises(ValueError, match='cursor'):
        widget.cursor = 'left'


@pytest.mark.parametrize('backend', ['bqplot', 'ginga'])
def test_sky_coordinate_format_property(backend):
    widget = _backend_widget(backend)
    assert widget.sky_coordinate_format == 'degrees'

    widget.sky_coordinate_format = 'sexagesimal'
    assert widget.sky_coordinate_format == 'sexagesimal'

    with pytest.raises(ValueError, match='sky_coordinate_format'):
        widget.sky_coordinate_format = 'radians'


@pytest.mark.parametrize('backend', ['bqplot', 'ginga'])
def test_mouse_move_updates_readout(backend, wcs):
    widget = _backend_widget(backend)
    data = np.arange(100.0 * 150).reshape(100, 150)
    widget.load_image(NDData(data=data, wcs=wcs))

    _simulate_mouse_move(widget, backend, 3.0, 2.0)
    expected = READOUT_TEMPLATE.format(
        format_cursor_text(3.0, 2.0, data=data, wcs=wcs))
    assert widget._cursor_readout.value == expected


@pytest.mark.parametrize('backend', ['bqplot', 'ginga'])
def test_mouse_move_before_load_does_not_raise(backend):
    widget = _backend_widget(backend)
    _simulate_mouse_move(widget, backend, 3.0, 2.0)
    assert widget._cursor_readout.value == 'Coordinates show up here'


@pytest.mark.parametrize('backend', ['bqplot', 'ginga'])
def test_sky_format_change_updates_readout_immediately(backend, wcs):
    # Regression test: changing sky_coordinate_format used to leave the
    # readout stale until the next mouse move.
    widget = _backend_widget(backend)
    data = np.arange(100.0 * 150).reshape(100, 150)
    widget.load_image(NDData(data=data, wcs=wcs))

    _simulate_mouse_move(widget, backend, 3.0, 2.0)
    widget.sky_coordinate_format = 'sexagesimal'
    expected = READOUT_TEMPLATE.format(
        format_cursor_text(3.0, 2.0, data=data, wcs=wcs,
                           sky_format='sexagesimal'))
    assert widget._cursor_readout.value == expected


@pytest.mark.parametrize('backend', ['bqplot', 'ginga'])
def test_sky_format_change_before_any_mouse_move(backend):
    widget = _backend_widget(backend)
    widget.sky_coordinate_format = 'sexagesimal'
    assert widget._cursor_readout.value == 'Coordinates show up here'


def test_backends_produce_identical_readout(wcs):
    # The payoff of the shared module: both backends must produce the
    # same readout string for the same cursor position and image.
    pytest.importorskip('bqplot')
    pytest.importorskip('ginga')

    data = np.arange(100.0 * 150).reshape(100, 150)
    readouts = []
    for backend in ('bqplot', 'ginga'):
        widget = _backend_widget(backend)
        widget.load_image(NDData(data=data, wcs=wcs))
        _simulate_mouse_move(widget, backend, 12.25, 34.75)
        readouts.append(widget._cursor_readout.value)

    assert readouts[0] == readouts[1]
    # floor(12.25 + 0.5) = 12
    assert f'X: {12:<5d}' in readouts[0]


def test_bqplot_cursor_assignment_regression():
    # Regression test: the AIDA slim-down (#215) removed the cursor
    # property from the bqplot backend, silently breaking
    # ``w.cursor = 'top'`` in the example notebook.
    pytest.importorskip('bqplot')
    widget = _backend_widget('bqplot')
    widget.cursor = 'top'
    assert widget.cursor == 'top'
