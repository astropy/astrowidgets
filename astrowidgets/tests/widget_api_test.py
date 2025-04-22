# TODO: How to enable switching out backend and still run the same tests?

import pytest

import numpy as np  # noqa: E402

from astropy.coordinates import SkyCoord  # noqa: E402
from astropy.io import fits  # noqa: E402
from astropy.nddata import NDData  # noqa: E402
from astropy.table import Table, vstack  # noqa: E402
from astropy import units as u  # noqa: E402
from astropy.wcs import WCS  # noqa: E402


class ImageWidgetAPITest:
    cursor_error_classes = (ValueError)

    @pytest.fixture
    def data(self):
        rng = np.random.default_rng(1234)
        return rng.random((100, 100))

    @pytest.fixture
    def wcs(self):
        # This is a copy/paste from the astropy 4.3.1 documentation...

        # Create a new WCS object.  The number of axes must be set
        # from the start
        w = WCS(naxis=2)

        # Set up an "Airy's zenithal" projection
        w.wcs.crpix = [-234.75, 8.3393]
        w.wcs.cdelt = np.array([-0.066667, 0.066667])
        w.wcs.crval = [0, -90]
        w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
        w.wcs.set_pv([(2, 1, 45.0)])
        return w

    # This setup is run before each test, ensuring that there are no
    # side effects of one test on another
    @pytest.fixture(autouse=True)
    def setup(self):
        """
        Subclasses MUST define ``image_widget_class`` -- doing so as a
        class variable does the trick.
        """
        self.image = self.image_widget_class(image_width=250, image_height=100)

    def test_width_height(self):
        assert self.image.image_width == 250
        assert self.image.image_height == 100

        width = 200
        height = 300
        self.image.image_width = width
        self.image.image_height = height
        assert self.image.image_width == width
        assert self.image.image_height == height

    def test_load_fits(self, data):
        hdu = fits.PrimaryHDU(data=data)
        self.image.load_fits(hdu)

    def test_load_nddata(self, data):
        nddata = NDData(data)
        self.image.load_nddata(nddata)

    def test_load_array(self, data):
        self.image.load_array(data)

    def test_center_on(self):
        self.image.center_on((10, 10))  # X, Y

    def test_offset_by(self, data, wcs):
        self.image.offset_by(10, 10)  # dX, dY

        # Testing offset by WCS requires a WCS. The viewer will (or ought to
        # have) taken care of setting up the WCS internally if initialized with
        # an NDData that has a WCS.
        ndd = NDData(data=data, wcs=wcs)
        self.image.load_nddata(ndd)

        self.image.offset_by(10 * u.arcmin, 10 * u.arcmin)

        # A mix of pixel and sky should produce an error
        with pytest.raises(ValueError, match='but dy is of type'):
            self.image.offset_by(10 * u.arcmin, 10)

        # A mix of inconsistent units should produce an error
        with pytest.raises(u.UnitConversionError):
            self.image.offset_by(1 * u.arcsec, 1 * u.AA)

    def test_zoom_level(self, data):
        # Set data first, since that is needed to determine zoom level
        self.image.load_array(data)
        self.image.zoom_level = 5
        assert self.image.zoom_level == 5

    def test_zoom(self):
        self.image.zoom_level = 3
        self.image.zoom(2)
        assert self.image.zoom_level == 6  # 3 x 2

    def test_marking_operations(self):
        marks = self.image.get_markers(marker_name="all")
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
            t = self.image.get_markers(marker_name='markymark')
        assert t is None
        assert len(warning_lines) == 1

        self.image.click_drag = True
        self.image.start_marking()
        assert not self.image.click_drag

        # Simulate a mouse click to add default marker name to the list.
        try:
            self.image._mouse_click_cb(self.image.viewer, None, 50, 50)
            assert self.image.get_marker_names() == [self.image._interactive_marker_set_name, 'markymark']
        except AttributeError:
            pass

        # Clear markers to not pollute other tests.
        self.image.stop_marking(clear_markers=True)

        assert self.image.is_marking is False
        assert self.image.get_markers(marker_name="all") is None

        # Hate this, should add to public API
        marknames = self.image._marktags
        assert len(marknames) == 0

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

        marknames = self.image._marktags
        assert marknames == set(['test1', 'test2'])
        # assert self.image.get_marker_names() == ['test1', 'test2']

        # No guarantee markers will come back in the same order, so sort them.
        t1 = self.image.get_markers(marker_name='test1')
        # Sort before comparing
        t1.sort('x')
        tab.sort('x')
        assert np.all(t1['x'] == tab['x'])
        assert (t1['y'] == tab['y']).all()

        # That should have given us two copies of the input table
        t2 = self.image.get_markers(marker_name="all")
        expected = vstack([tab, tab], join_type='exact')
        # Sort before comparing
        t2.sort(['x', 'y'])
        expected.sort(['x', 'y'])
        assert (t2['x'] == expected['x']).all()
        assert (t2['y'] == expected['y']).all()

        self.image.remove_markers(marker_name='test1')
        marknames = self.image._marktags
        assert marknames == set(['test2'])
        # assert self.image.get_marker_names() == ['test2']

        # Ensure unable to mark with reserved name
        for name in self.image.RESERVED_MARKER_SET_NAMES:
            with pytest.raises(ValueError, match='not allowed'):
                self.image.add_markers(tab, marker_name=name)

        # Add markers with no marker name and check we can retrieve them
        # using the default marker name
        self.image.add_markers(tab, x_colname='x', y_colname='y',
                               skycoord_colname='coord')
        # Don't care about the order of the marker names so use set instead of
        # list.
        marknames = self.image._marktags
        assert (set(marknames) == set(['test2', self.image._default_mark_tag_name]))
        # assert (set(self.image.get_marker_names()) ==
        #         set(['test2', self.image._default_mark_tag_name]))

        # Clear markers to not pollute other tests.
        self.image.reset_markers()
        marknames = self.image._marktags
        assert len(marknames) == 0
        assert self.image.get_markers(marker_name="all") is None
        with pytest.warns(UserWarning, match='is empty'):
            assert self.image.get_markers(marker_name=self.image._default_mark_tag_name) is None

        with pytest.raises(ValueError, match="No markers named 'test1'"):
            self.image.get_markers(marker_name='test1')
        with pytest.raises(ValueError, match="No markers named 'test2'"):
            self.image.get_markers(marker_name='test2')

    def test_get_markers_accepts_list_of_names(self):
        # Check that the get_markers method accepts a list of marker names
        # and returns a table with all the markers from all the named sets.
        rng = np.random.default_rng(1234)
        data = rng.integers(0, 100, (5, 2))
        tab = Table(data=data, names=['x', 'y'])
        self.image.add_markers(tab, x_colname='x', y_colname='y',
                                 skycoord_colname='coord', marker_name='test1')
        self.image.add_markers(tab, x_colname='x', y_colname='y',
                                    skycoord_colname='coord', marker_name='test2')

        # No guarantee markers will come back in the same order, so sort them.
        t1 = self.image.get_markers(marker_name=['test1', 'test2'])
        # Sort before comparing
        t1.sort('x')
        expected = vstack([tab, tab], join_type='exact')
        expected.sort('x')
        assert np.all(t1['x'] == expected['x'])
        assert (t1['y'] == expected['y']).all()

    def test_remove_markers(self):
        with pytest.raises(ValueError, match='arf'):
            self.image.remove_markers(marker_name='arf')

    def test_remove_markers_name_all(self):
        rng = np.random.default_rng(1234)
        data = rng.integers(0, 100, (5, 2))
        tab = Table(data=data, names=['x', 'y'])
        self.image.add_markers(tab, marker_name='test1')
        self.image.add_markers(tab, marker_name='test2')

        self.image.remove_markers(marker_name='all')
        assert self.image.get_markers(marker_name='all') is None

    def test_adding_markers_as_world(self, data, wcs):
        ndd = NDData(data=data, wcs=wcs)
        self.image.load_nddata(ndd)

        # Add markers using world coordinates
        rng = np.random.default_rng(9435)

        pixels = rng.integers(0, 100, (5, 2))
        marks_pix = Table(data=pixels, names=['x', 'y'], dtype=('float', 'float'))
        marks_world = wcs.pixel_to_world(marks_pix['x'], marks_pix['y'])
        marks_coords = SkyCoord(marks_world, unit='degree')
        mark_coord_table = Table(data=[marks_coords], names=['coord'])
        self.image.add_markers(mark_coord_table, use_skycoord=True)
        result = self.image.get_markers()
        # Check the x, y positions as long as we are testing things...
        np.testing.assert_allclose(result['x'], marks_pix['x'])
        np.testing.assert_allclose(result['y'], marks_pix['y'])
        np.testing.assert_allclose(result['coord'].ra.deg,
                                   mark_coord_table['coord'].ra.deg)
        np.testing.assert_allclose(result['coord'].dec.deg,
                                   mark_coord_table['coord'].dec.deg)

    def test_stretch(self):
        original_stretch = self.image.stretch

        with pytest.raises(ValueError, match='must be one of'):
            self.image.stretch = 'not a valid value'

        # A bad value should leave the stretch unchanged
        assert self.image.stretch is original_stretch

        self.image.stretch = 'log'
        # A valid value should change the stretch
        assert self.image.stretch is not original_stretch

    def test_cuts(self, data):
        with pytest.raises(ValueError, match='must be one of'):
            self.image.cuts = 'not a valid value'

        with pytest.raises(ValueError, match='must have length 2'):
            self.image.cuts = (1, 10, 100)

        assert 'histogram' in self.image.autocut_options

        # Setting using histogram requires data
        self.image.load_array(data)
        self.image.cuts = 'histogram'
        assert len(self.image.cuts) == 2

        self.image.cuts = (10, 100)
        assert self.image.cuts == (10, 100)

    def test_colormap(self):
        cmap_desired = 'gray'
        cmap_list = self.image.colormap_options
        assert len(cmap_list) > 0 and cmap_desired in cmap_list
        self.image.set_colormap(cmap_desired)

    def test_cursor(self):
        assert self.image.cursor in self.image.ALLOWED_CURSOR_LOCATIONS
        with pytest.raises(self.cursor_error_classes):
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
        self.image.click_drag = False
        self.image._is_marking = True
        with pytest.raises(ValueError, match=r'([Ii]nteractive marking)|(while in marking mode)'):
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
        with pytest.raises(ValueError, match=r'([Ii]nteractive marking)|(while in marking mode)'):
            self.image.click_center = True
        self.image._is_marking = False

    def test_scroll_pan(self):
        # Make sure scroll_pan is actually settable
        for value in [True, False]:
            self.image.scroll_pan = value
            assert self.image.scroll_pan is value

    def test_save(self, tmp_path):
        filename = tmp_path / 'woot.png'
        self.image.save(filename)

    def test_save_overwrite(self, tmp_path):
        filename = tmp_path / 'woot.png'

        # First write should be fine
        self.image.save(filename)

        # Second write should raise an error because file exists
        with pytest.raises(FileExistsError):
            self.image.save(filename)

        # Using overwrite should save successfully
        self.image.save(filename, overwrite=True)
