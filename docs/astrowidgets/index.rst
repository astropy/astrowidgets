Image widget for Jupyter Lab/notebook
=====================================

Getting started
---------------

Make a viewer
+++++++++++++

The snippet below is all you need to make an image widget. The widget is part
of the `ipywidgets framework <https://ipywidgets.rtfd.io>`_ so that it can
be easily integrated with other controls::

    >>> from astrowidgets import ImageWidget
    >>> image = ImageWidget()
    >>> display(image)

Loading an image
++++++++++++++++

An empty viewer is not very useful, though, so load some data from a FITS
file. The FITS file at the link below is an image of the field of the exoplanet
Kelt-16, and also contains part of the Veil Nebula::

    >>> image.load_fits('https://zenodo.org/record/3356833/files/kelt-16-b-S001-R001-C084-r.fit.bz2?download=1')

The image widget can also load a Numpy array via
`~astrowidgets.ImageWidget.load_array`. It also understands astropy
`~astropy.nddata.NDData` objects; load them via
`~astrowidgets.ImageWidget.load_data`.

Navigation
++++++++++

In the default configuration, basic navigation is done using these controls:

* scroll to pan
* use ``+``/``-`` to zoom in/out (cursor must be over the image for this to work)
* right-click and drag to change contrast DS9-style

API
+++

One important design goal is to make all functionality available by a compact,
clear API. The `target API <https://github.com/eteq/nb-astroimage-api>`_ still
needs a few features (e.g., blink), but much of it is already implemented.

The API-first approach means that manipulating the view programmatically is straightforward.
For example, centering on the position of the object, Kelt-16, and zooming in to 8x the natural
pixel scale is straightforward::

    >>> from astropy.coordinates import SkyCoord
    >>> image.center_on(SkyCoord.from_name('kelt-16'))
    >>> image.zoom_level = 8

A more detailed description of the interface and the :ref:`api-docs` are available.

.. toctree::
    :maxdepth: 2

    api.rst

Example Notebooks
-----------------

* `astrowidgets using the Ginga backend <https://github.com/astropy/astrowidgets/blob/main/example_notebooks/ginga_widget.ipynb>`_
* `Using named markers to keep track of logically related markers <https://github.com/astropy/astrowidgets/blob/main/example_notebooks/named_markers.ipynb>`_
* `Demonstration of GUI interactions <https://github.com/astropy/astrowidgets/blob/main/example_notebooks/gui_interactions.ipynb>`_


