# This file is used to configure the behavior of pytest when using the Astropy
# test infrastructure.
try:
    from pytest_astropy_header.display import (PYTEST_HEADER_MODULES,
                                               TESTED_VERSIONS)
except ImportError:
    PYTEST_HEADER_MODULES = {}
    TESTED_VERSIONS = {}

# Uncomment and customize the following lines to add/remove entries from
# the list of packages for which version numbers are displayed when running
# the tests. Making it pass for KeyError is essential in some cases when
# the package uses other astropy affiliated packages.
PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
PYTEST_HEADER_MODULES['Ginga'] = 'ginga'
PYTEST_HEADER_MODULES.pop('h5py', None)
PYTEST_HEADER_MODULES.pop('Pandas', None)

# Uncomment the following lines to display the version number of the
# package rather than the version number of Astropy in the top line when
# running the tests.
try:
    from astrowidgets import __version__ as version
except ImportError:
    version = 'unknown'
TESTED_VERSIONS['astrowidgets'] = version
