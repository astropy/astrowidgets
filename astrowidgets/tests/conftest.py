import numpy as np
import pytest

from astropy.wcs import WCS


@pytest.fixture
def wcs():
    # The Airy-projection example WCS from the astropy documentation, the
    # same one used by the wcs fixture in
    # astro_image_display_api.api_test.ImageAPITest.
    w = WCS(naxis=2)
    w.wcs.crpix = [-234.75, 8.3393]
    w.wcs.cdelt = np.array([-0.066667, 0.066667])
    w.wcs.crval = [0, -90]
    w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    w.wcs.set_pv([(2, 1, 45.0)])
    return w
