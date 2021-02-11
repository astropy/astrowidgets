# TODO: How to enable switching out backend and still run the same tests?

import pytest

ginga = pytest.importorskip("ginga")

import numpy as np  # noqa: E402

from astropy.io import fits  # noqa: E402
from astropy.nddata import NDData  # noqa: E402
from astropy.table import Table, vstack  # noqa: E402

from ginga.ColorDist import ColorDistBase  # noqa: E402
from ginga.misc.log import get_logger  # noqa: E402

from astrowidgets.ginga import ImageWidget  # noqa: E402


class TestGingaWidget:
    def setup_class(self):
        rng = np.random.default_rng(1234)
        self.data = rng.random((100, 100))

        logger = get_logger('my_viewer', log_stderr=False, level=30)
        self.image = ImageWidget(logger, image_width=250, image_height=100)

    def test_width_height(self):
        assert self.image.image_width == 250
        assert self.image.image_height == 100

        width = 200
        height = 300
        self.image.image_width = width
        self.image.image_height = height
        assert self.image.image_width == width
        assert self.image.image_height == height
        assert self.image.viewer.get_window_size() == (width, height)

    def test_load_fits(self):
        hdu = fits.PrimaryHDU(data=self.data)
        self.image.load_fits(hdu)

    def test_load_nddata(self):
        nddata = NDData(self.data)
        self.image.load_nddata(nddata)

    def test_load_array(self):
        self.image.load_array(self.data)

    def test_center_on(self):
        self.image.center_on((10, 10))  # X, Y

    def test_offset_to(self):
        self.image.offset_to(10, 10)  # dX, dY

    def test_zoom_level(self):
        self.image.zoom_level = 5
        assert self.image.zoom_level == 5

    def test_zoom(self):
        self.image.zoom_level = 3
        self.image.zoom(2)
        assert self.image.zoom_level == 6  # 3 x 2

    @pytest.mark.xfail(reason='Not implemented yet')
    def test_select_points(self):
        self.image.select_points()

    def test_marking_operations(self):
        marks = self.image.get_all_markers()
        assert marks is None
        assert not self.image.is_marking

        # Ensure you cannot set it like this.
        with pytest.raises(AttributeError):
            self.image.is_marking = True

        # Setting these to check that start_marking affects them.
        self.image.click_center = True  # Disables click_drag
        assert self.image.click_center
        self.image.scroll_pan = False
        assert not self.image.scroll_pan

        # Set the marker style
        marker_style = {'color': 'yellow', 'radius': 10, 'type': 'cross'}
        m_str = str(self.image.marker)
        for key in marker_style.keys():
            assert key in m_str

        self.image.start_marking(marker_name='markymark', marker=marker_style)
        assert self.image.is_marking
        assert self.image.marker == marker_style
        assert not self.image.click_center
        assert not self.image.click_drag

        # scroll_pan better activate when marking otherwise there is
        # no way to pan while interactively marking
        assert self.image.scroll_pan

        # Make sure that when we stop_marking we get our old controls back.
        self.image.stop_marking()
        assert self.image.click_center
        assert not self.image.click_drag
        assert not self.image.scroll_pan

        # Regression test for GitHub Issue 97:
        # Marker name with no markers should give warning.
        with pytest.warns(UserWarning, match='is empty') as warning_lines:
            t = self.image.get_markers_by_name('markymark')
        assert t is None
        assert len(warning_lines) == 1

        self.image.click_drag = True
        self.image.start_marking()
        assert not self.image.click_drag

        # Simulate a mouse click to add default marker name to the list.
        self.image._mouse_click_cb(self.image.viewer, None, 50, 50)
        assert self.image.get_marker_names() == [self.image._interactive_marker_set_name, 'markymark']

        # Clear markers to not pollute other tests.
        self.image.stop_marking(clear_markers=True)

        assert self.image.is_marking is False
        assert self.image.get_all_markers() is None
        assert len(self.image.get_marker_names()) == 0

        # Make sure that click_drag is restored as expected
        assert self.image.click_drag

    def test_add_markers(self):
        rng = np.random.default_rng(1234)
        data = rng.integers(0, 100, (5, 2))
        orig_tab = Table(data=data, names=['x', 'y'], dtype=('float', 'float'))
        tab = Table(data=data, names=['x', 'y'], dtype=('float', 'float'))
        self.image.add_markers(tab, x_colname='x', y_colname='y',
                               skycoord_colname='coord', marker_name='test1')

        # Make sure setting didn't change the default name
        assert self.image._default_mark_tag_name == 'default-marker-name'

        # Regression test for GitHub Issue 45:
        # Adding markers should not modify the input data table.
        assert (tab == orig_tab).all()

        # Add more markers under different name.
        self.image.add_markers(tab, x_colname='x', y_colname='y',
                               skycoord_colname='coord', marker_name='test2')
        assert self.image.get_marker_names() == ['test1', 'test2']

        # No guarantee markers will come back in the same order, so sort them.
        t1 = self.image.get_markers_by_name('test1')
        # Sort before comparing
        t1.sort('x')
        tab.sort('x')
        assert np.all(t1['x'] == tab['x'])
        assert (t1['y'] == tab['y']).all()

        # That should have given us two copies of the input table
        t2 = self.image.get_all_markers()
        expected = vstack([tab, tab], join_type='exact')
        # Sort before comparing
        t2.sort(['x', 'y'])
        expected.sort(['x', 'y'])
        assert (t2['x'] == expected['x']).all()
        assert (t2['y'] == expected['y']).all()

        self.image.remove_markers_by_name('test1')
        assert self.image.get_marker_names() == ['test2']

        # Ensure unable to mark with reserved name
        for name in self.image.RESERVED_MARKER_SET_NAMES:
            with pytest.raises(ValueError, match='not allowed'):
                self.image.add_markers(tab, marker_name=name)

        # Clear markers to not pollute other tests.
        self.image.remove_all_markers()
        assert len(self.image.get_marker_names()) == 0
        assert self.image.get_all_markers() is None
        assert self.image.get_markers_by_name(self.image._default_mark_tag_name) is None

        with pytest.raises(ValueError, match="No markers named 'test1'"):
            self.image.get_markers_by_name('test1')
        with pytest.raises(ValueError, match="No markers named 'test2'"):
            self.image.get_markers_by_name('test2')

    def test_remove_markers(self):
        with pytest.raises(ValueError, match='arf'):
            self.image.remove_markers_by_name('arf')

    def test_stretch(self):
        with pytest.raises(ValueError, match='must be one of'):
            self.image.stretch = 'not a valid value'

        self.image.stretch = 'log'
        assert isinstance(self.image.stretch, ColorDistBase)

    def test_cuts(self):
        with pytest.raises(ValueError, match='must be one of'):
            self.image.cuts = 'not a valid value'

        with pytest.raises(ValueError, match=r'must be given as \(low, high\)'):
            self.image.cuts = (1, 10, 100)

        assert 'histogram' in self.image.autocut_options

        self.image.cuts = 'histogram'
        np.testing.assert_allclose(
            self.image.cuts, (3.948844e-04, 9.990224e-01), rtol=1e-6)

        self.image.cuts = (10, 100)
        assert self.image.cuts == (10, 100)

    def test_colormap(self):
        cmap_desired = 'gray'
        cmap_list = self.image.colormap_options
        assert len(cmap_list) > 0 and cmap_desired in cmap_list
        self.image.set_colormap(cmap_desired)

    def test_cursor(self):
        assert self.image.cursor in self.image.ALLOWED_CURSOR_LOCATIONS
        with pytest.raises(ValueError):
            self.image.cursor = 'not a valid option'
        self.image.cursor = 'bottom'
        assert self.image.cursor == 'bottom'

    def test_click_drag(self):
        # Set this to ensure that click_drag turns it off
        self.image.click_center = True

        # Make sure that setting click_drag to False does not turn off
        # click_center.
        self.image.click_drag = False
        assert self.image.click_center

        self.image.click_drag = True
        assert not self.image.click_center

        # If is_marking is true then trying to enable click_drag should fail
        self.image._is_marking = True
        self.image.click_drag = False
        with pytest.raises(ValueError, match='Interactive marking'):
            self.image.click_drag = True
        self.image._is_marking = False

    def test_click_center(self):
        # Set this to ensure that click_center turns it off
        self.image.click_drag = True

        # Make sure that setting click_center to False does not turn off
        # click_draf.
        self.image.click_center = False
        assert self.image.click_drag

        self.image.click_center = True
        assert not self.image.click_drag

        # If is_marking is true then trying to enable click_center should fail
        self.image._is_marking = True
        self.image.click_center = False
        with pytest.raises(ValueError, match='Interactive marking'):
            self.image.click_center = True
        self.image._is_marking = False

    def test_scroll_pan(self):
        # Make sure scroll_pan is actually settable
        for value in [True, False]:
            self.image.scroll_pan = value
            assert self.image.scroll_pan is value

    def test_save(self, tmpdir):
        with pytest.raises(ValueError, match='not supported'):
            self.image.save(str(tmpdir.join('woot.jpg')))

        filename = str(tmpdir.join('woot.png'))
        self.image.save(filename)

        with pytest.raises(ValueError, match='exists'):
            self.image.save(filename)

        self.image.save(filename, overwrite=True)
