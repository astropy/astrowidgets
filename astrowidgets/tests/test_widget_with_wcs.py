# TODO: How to enable switching out backend and still run the same tests?

import pytest

ginga = pytest.importorskip("ginga")

import numpy as np  # noqa: E402
from astropy.coordinates import SkyCoord  # noqa: E402
from astropy.nddata import CCDData  # noqa: E402
from astropy.table import Table  # noqa: E402
from astropy.wcs import WCS  # noqa: E402

from ginga.misc.log import get_logger  # noqa: E402

from astrowidgets.ginga import ImageWidget  # noqa: E402


def _make_fake_ccd(with_wcs=True):
    """Generate a CCDData object for use with ImageWidget tests.

    Parameters
    ----------
    with_wcs : bool, optional
        If `True` the image will have a WCS attached to it,
        which is useful for some of the marker tests.

    Returns
    -------
    image : `astropy.nddata.CCDData`
        CCD image.

    """
    rng = np.random.default_rng(1234)
    npix_side = 100
    fake_image = rng.random((npix_side, npix_side))
    if with_wcs:
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = (fake_image.shape[0] / 2, fake_image.shape[1] / 2)
        wcs.wcs.ctype = ('RA---TAN', 'DEC--TAN')
        wcs.wcs.crval = (314.275419158, 31.6662781301)
        wcs.wcs.pc = [[0.000153051015113, -3.20700931602e-05],
                      [3.20704370872e-05, 0.000153072382405]]
    else:
        wcs = None

    return CCDData(data=fake_image, wcs=wcs, unit='adu')


class TestGingaWidgetWithWCS:
    def setup_class(self):
        self.fake_ccd = _make_fake_ccd(with_wcs=True)

        logger = get_logger('my_viewer', log_stderr=False, level=30)
        self.image = ImageWidget(logger, pixel_coords_offset=0)
        self.image.load_nddata(self.fake_ccd)

    def test_adding_markers_as_world_recovers_with_get_markers(self):
        """Make sure that our internal conversion from world to pixel
        coordinates doesn't mess anything up.
        """
        npix_side = self.fake_ccd.shape[0]
        wcs = self.fake_ccd.wcs

        # Get me 100 positions please, not right at the edge
        rng = np.random.default_rng(1234)
        marker_locs = rng.integers(10, high=npix_side - 10, size=(100, 2))
        marks_pix = Table(data=marker_locs, names=('x', 'y'))
        marks_world = wcs.all_pix2world(marker_locs, 0)
        marks_coords = SkyCoord(marks_world, unit='degree')
        mark_coord_table = Table(data=[marks_coords], names=('coord', ))
        self.image.add_markers(mark_coord_table, use_skycoord=True)

        result = self.image.get_all_markers()
        np.testing.assert_allclose(
            result['coord'].ra.deg, mark_coord_table['coord'].ra.deg)
        np.testing.assert_allclose(
            result['coord'].dec.deg, mark_coord_table['coord'].dec.deg)
        # Might as well check the x, y positions too.
        np.testing.assert_allclose(result['x'], marks_pix['x'])
        np.testing.assert_allclose(result['y'], marks_pix['y'])

        # Clear markers before running other tests.
        self.image.remove_all_markers()

    def test_empty_marker_name_works_with_all(self):
        """Regression test for GitHub Issue 97.

        This particular test checks that an empty table doesn't break
        get_all_markers(). The bug only comes up if there is a coordinate
        column, so use a fake image a WCS.

        """
        x = np.array([20, 30, 40])
        y = np.array([40, 80, 100])
        input_markers = Table(data=[x, y], names=['x', 'y'])
        # Add some markers with our own name
        self.image.add_markers(input_markers, marker_name='nowcs')

        # Start marking to create a new marker set that is empty
        self.image.start_marking(marker_name='empty')
        self.image.stop_marking()

        assert self.image.get_marker_names() == ['empty', 'nowcs']

        with pytest.warns(UserWarning, match='is empty'):
            marks = self.image.get_all_markers()
        assert len(marks) == len(x)
        assert 'empty' not in marks['marker name']

        # Clear markers before running other tests.
        self.image.remove_all_markers()
