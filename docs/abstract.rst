.. _abstract_widget_intro:

Understanding BaseImageWidget
=============================

``astrowidgets`` provides an abstract class called
`~astrowidgets.core.BaseImageWidget` to allow developers from different
visualization tools (hereafter known as "backends") to implement their own
solutions using the same set of API. This design is based on
`nb-astroimage-api <https://github.com/eteq/nb-astroimage-api>`_, with the
goal of making all functionality available by a compact and clear API.
This API-first approach would allow manipulating the view programmatically
in a reproducible way.

The idea of the abstract class is that ``astrowidgets`` users would be able
to switch to the backend of their choice without much refactoring of their own.
This would enable, say, astronomers with different backend preferences to
collaborate more easily via Jupyter Lab/Notebook.

The following sub-sections lay out the envisioned high-level usage of
``astrowidgets`` regardless of backend inside Jupyter Lab/Notebook.
However, the examples are not exhaustive. For the full API definition,
please see :ref:`abstract_api`.

.. _abstract_viewer:

Creating a Viewer
-----------------

The snippet below is all you should need to make an image widget.
The widget should be a part of the
`ipywidgets framework <https://ipywidgets.rtfd.io>`_ so that it can
be easily integrated with other controls:

.. code-block:: python

    from astrowidgets.somebackend import ImageWidget
    image = ImageWidget()
    image

.. _abstract_image_load:

Loading an Image
----------------

To load data into the empty viewer created in :ref:`abstract_viewer`,
there should be methods to load different formats:

.. code-block:: python

    # FITS image of the field of the exoplanet Kelt-16,
    # and also contains part of the Veil Nebula
    filename = 'https://zenodo.org/record/3356833/files/kelt-16-b-S001-R001-C084-r.fit.bz2?download=1'
    image.load_fits(filename)

.. code-block:: python

    # A Numpy array
    import numpy as np
    arr = np.arange(100).reshape(10, 10)
    image.load_array(arr)

.. code-block:: python

    # An astropy.nddata.NDData object
    from astropy.io import fits
    from astropy.nddata import NDData
    from astropy.wcs import WCS
    with fits.open(filename) as pf:
        data = NDData(pf[0].data, wcs=WCS(pf[0].header))
    image.load_nddata(data)

If additional format support is desired, the API could be added to the
abstract base class if the new format is widely supported and not specific
to a certain backend implementation.

.. _abstract_cursor_info:

Cursor Info Display
-------------------

The widget actually consists of two child widgets:

* The image display.
* Cursor information panel with the following:
    * X and Y cursor locations, taking
      `~astrowidgets.core.BaseImageWidget.pixel_offset` into account.
    * RA and Dec calculated from the cursor location using the image's WCS,
      if available. It is up to the backend on how to handle WCS projection
      or distortion.
    * Value of the image under the cursor.

The cursor information panel can have three different states:

* Positioned below the image display.
* Positioned above the image display.
* Not displayed.

This state can be set using the `~astrowidgets.core.BaseImageWidget.cursor`
property.

.. _abstract_size:

Changing Display Size
---------------------

There should be a programmatic way to change the display size of the display
widget:

.. code-block:: python

    # The height would auto-adjust
    image.image_width = 500  # pixels

    # The width would auto-adjust
    image.image_height = 500  # pixels

.. _abstract_colormap:

Changing Colormap
-----------------

There should be a programmatic way to change the colormap of the display.
However, the available colormaps may differ from backend to backend:

.. code-block:: python

    image.set_colormap('viridis')

.. _abstract_controls:

Mouse/Keyboard Controls
-----------------------

Mouse interaction using clicks and scroll should be supported.
Keyboard controls would also be desirable. These controls should be active
when cursor is over the display, but not otherwise.
For example, but not limited to:

* Scrolling to pan up/down the image.
* Using ``+``/``-`` to zoom in/out.
* Using click-and-drag to change the contrast of the image.

In the event where the same click/button can be overloaded, the active
functionality can be controlled by the following properties:

* `~astrowidgets.core.BaseImageWidget.click_center`
* `~astrowidgets.core.BaseImageWidget.click_drag`
* `~astrowidgets.core.BaseImageWidget.scroll_pan`

There should be programmatic ways to perform these controls as well:

.. code-block:: python

    # Centering on sky coordinates
    from astropy.coordinates import SkyCoord
    image.center_on(SkyCoord.from_name('kelt-16'))

    # Centering on pixel coordinates
    image.center_on((100, 100))

    # Moving the center using sky coordinates
    from astropy import units as u
    image.offset_to(0.1 * u.arcsec, 0.1 * u.arcsec, skycoord_offset=True)

    # Moving the center by pixels
    image.offset_to(10, 10)

    # Zooming (two different ways)
    image.zoom(2)
    image.zoom_level = 1

    # Changing the display stretch
    image.stretch = 'log'

    # Changing the cut levels (two different ways)
    image.cuts = 'histogram'
    image.cuts = (0, 10)  # (low, high)

Please also see :ref:`abstract_marking`.

.. _abstract_marking:

Marking Objects
---------------

Another important aspect is to allow users to either interactively or
programmatically mark objects of interest on the displayed image.
Marking mode is tracked using the
`~astrowidgets.core.BaseImageWidget.is_marking`
property and can be turned on and off using
:meth:`~astrowidgets.core.BaseImageWidget.start_marking` and
:meth:`~astrowidgets.core.BaseImageWidget.stop_marking`, respectively.
The marker appearance can be changed using
`~astrowidgets.core.BaseImageWidget.marker`.

For interactive marking, after a user runs ``start_marking`` but before
``stop_marking``, a click on the image display would mark the object under
the cursor.

For programmatic marking, user can first build a `~astopy.table.Table` with
either pixel or sky coordinates, and then pass it into
:meth:`~astrowidgets.core.BaseImageWidget.add_markers`.

User can then call
:meth:`~astrowidgets.core.BaseImageWidget.get_markers_by_name` or
:meth:`~astrowidgets.core.BaseImageWidget.get_all_markers` to obtain the
marked locations.

To remove the markers, user can call
:meth:`~astrowidgets.core.BaseImageWidget.remove_markers_by_name` or
:meth:`~astrowidgets.core.BaseImageWidget.remove_all_markers`, as appropriate.

To put this all together, here is an example workflow (out of many)
that may happen:

1. User calls ``start_marking`` to begin the interactive marking session.
2. User clicks on two stars.
3. User calls ``stop_marking`` to end the interactive marking session.
4. User reads a table from a collaborator containing several galaxies in the
   field of view.
5. User changes the marker style from a red circle to a green square by
   modifying the ``marker`` property.
6. User programmatically marks the galaxies on display with the new marker style
   and a new marker name using ``add_markers``.
7. User obtains all the marked locations for post-processing using
   ``get_all_markers``.
8. User removes all the markers from display using ``remove_all_markers``.

.. _abstract_save:

Saving an Image
---------------

The image display can be programmatically saved to a file, but not the
:ref:`abstract_cursor_info`. Supported output format is controlled by the
individual backend. For example:

.. code-block:: python

    image.save('myimage.png')

.. _example_notebooks:

Example Notebooks
-----------------

Please see the `example notebooks folder <https://github.com/astropy/astrowidgets/blob/master/example_notebooks/>`_
for examples using a concrete implementation of this abstract class.
Backend-dependent dependencies are required to run them.

.. _abstract_api:

API
---

.. automodapi:: astrowidgets
    :no-inheritance-diagram:
