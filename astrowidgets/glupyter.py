"""Module containing ``astrowidgets`` implemented with ``glue-jupyter``
(that uses ``bqplot``) backend.

For this to work, ``astrowidgets`` must be installed along with the optional
dependencies specified for the ``glue-jupyter`` (a.k.a. Glupyter or
``glupyter``) backend; e.g.,::

    pip install 'astrowidgets[glupyter]'

"""
from glue_jupyter import jglue

from astrowidgets.core import BaseImageWidget

__all__ = ['ImageWidget']


class ImageWidget(BaseImageWidget):
    """Image widget for Jupyter notebook using ``glue-jupyter``/``bqplot``
    viewer.

    Parameters
    ----------
    kwargs : dict
        See `~astrowidgets.core.BaseImageWidget`.

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self._app = jglue()
        self._viewer = self._app.imshow()

        # UNTIL HERE -- imshow errors out without data
        # TODO: This builds on PR 126. Example notebook in test_data/ztmp...

    @property
    def viewer(self):
        self._viewer
