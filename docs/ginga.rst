.. _ginga_backend:

Widget with Ginga Toolkit
=========================

``astrowidgets`` comes with an example concrete implementation using
`Ginga <https://ginga.readthedocs.io>`_ as a backend:

.. code-block:: python

    from astrowidgets.ginga import ImageWidget
    from ginga.misc.log import get_logger
    logger = get_logger('my_viewer', log_stderr=False, log_file='ginga.log',
                        level=40)
    image = ImageWidget(logger)

Please see the `Ginga example notebooks folder <https://github.com/astropy/astrowidgets/blob/master/example_notebooks/ginga/>`_
for examples using this implementation.

.. _ginga_dependencies:

Dependencies
------------

The following dependecies need to be installed separately if you wish to use
the Ginga implementation:

* ``ginga>=2.7.1``
* ``pillow``
* ``freetype``
* ``aggdraw``
* ``opencv`` (optional, not required but will improve performance)

.. note::

    For vectorized drawing in ``aggdraw``, you can clone
    https://github.com/ejeschke/aggdraw/ and install its ``vectorized-drawing``
    branch from source.

For Windows Users
^^^^^^^^^^^^^^^^^

It is a known issue that ``FREETYPE_ROOT`` is not set properly if you do
``conda install aggdraw`` on Windows
(https://github.com/conda-forge/freetype-feedstock/issues/12), which results
in ``aggdraw cannot load font (no text renderer)`` error message when
using the widget with Ginga toolkit. The solution is to update to ``aggdraw``
1.3.5 or later; e.g., ``conda install aggdraw=1.3.5``.

.. _ginga_opencv:

Using OpenCV
------------

If you wish to use `OpenCV <https://docs.opencv.org/master/index.html>`_
to handle the drawing in Ginga, you have two options:

Install OpenCV with pip
^^^^^^^^^^^^^^^^^^^^^^^

If you are using pip it looks like the best option is to use the
``opencv-python`` package, which provides pre-built binaries of most of OpenCV::

    pip install opencv-python

However, the `opencv-python project <https://github.com/skvark/opencv-python>`_
is quite clear about being "unofficial" so you should probably read about
the project before using.

Install OpenCV with conda
^^^^^^^^^^^^^^^^^^^^^^^^^

This should work on conda::

    conda install -c conda-forge opencv

If, after installing ``opencv``, you get a warning like this::

    UserWarning: install opencv or set use_opencv=False

Then, you should try installing a newer version of ``freetype``::

    conda install 'freetype\>=2.10'

For more details, see `this discussion of opencv and astrowidgets
<https://github.com/astropy/astrowidgets/issues/90>`_.

.. _ginga_imagewidget_api:

API
---

.. automodapi:: astrowidgets.ginga
