Installation
============

Install with pip
----------------

Note that to use ``astrowidgets`` with Jupyter Lab you need to
install `nodejs from here <https://nodejs.org/en/download/>`_::

    pip install astrowidgets
    #
    # ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
    # INSTALL NODEJS FROM https://nodejs.org/en/download/
    # ↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
    #
    # (the nodejs on PyPI is *NOT* the right thing)
    #
    jupyter labextension install @jupyter-widgets/jupyterlab-manager

Install with conda
------------------

conda installation::

    conda install -c conda-forge astrowidgets nodejs
    jupyter labextension install @jupyter-widgets/jupyterlab-manager
    # Optional, but avoids a rebuild when you open lab
    jupyter lab build

.. note::

    If you are using a virtual environment and
    are experiencing problems with the installation, be sure to check and use
    the suitable ``pip`` version. Also, if you are having problem even with the following
    dependencies installed, please make sure you have the most up-to-date versions.
    Upgrade dependencies when and where necessary.

The following dependencies are needed *in the kernel of execution*
to use ``astrowidgets`` in either Jupyter Lab or Notebook. They should be installed
automatically when you install astrowidgets:

* ``python>=3.8``
* ``numpy``
* ``astropy``
* ``ipywidgets>=7.5``
* ``ipyevents>=0.6.3``
* ``ginga>=2.7.1``
* ``pillow``
* ``freetype``
* ``aggdraw``
* ``jupyterlab>=3``
* ``nodejs``
* ``opencv`` (optional, not installed by default)

After installing dependencies, for Jupyter Lab, run::

    jupyter labextension install @jupyter-widgets/jupyterlab-manager

For those using ``conda``, dependencies from the ``conda-forge`` channel
should be sufficient unless stated otherwise.

Using OpenCV
------------

If you wish to use `OpenCV <https://docs.opencv.org/master/index.html>`_ to handle the
drawing in Ginga, you have two options:

Install OpenCV with pip
^^^^^^^^^^^^^^^^^^^^^^^

If you are using pip it looks like the best option is to use the
``opencv-python`` package, which provides pre-built binaries of most of OpenCV::

    pip install opencv-python

However, the `opencv-python project
<https://github.com/skvark/opencv-python>`_ is quite clear about being
"unofficial" so you should probably read about the project before using.

Install OpenCV with conda
^^^^^^^^^^^^^^^^^^^^^^^^^

This should work on conda::

    conda install -c conda-forge opencv

If, after installing ``opencv``, you get a warning like this::

    astrowidgets/core.py:72: UserWarning: install opencv or set use_opencv=False
    warnings.warn('install opencv or set use_opencv=False')

then you should try installing a newer version of ``freetype``::

    conda install 'freetype\>=2.10'

For more details, see `this discussion of opencv and astrowidgets
<https://github.com/astropy/astrowidgets/issues/90>`_.

Widget with Ginga toolkit
-------------------------

.. note::

    For vectorized drawing in ``aggdraw``, you can clone
    https://github.com/ejeschke/aggdraw/ and install its ``vectorized-drawing``
    branch from source.


Notes for Windows users
-----------------------

aggdraw
^^^^^^^

It is a known issue that ``FREETYPE_ROOT`` is not set properly if you do
``conda install aggdraw`` on Windows
(https://github.com/conda-forge/freetype-feedstock/issues/12), which results
in ``aggdraw cannot load font (no text renderer)`` error message when
using the widget with Ginga toolkit. The solution is to update to ``aggdraw``
1.3.5 or later; e.g., ``conda install aggdraw=1.3.5``.

nodejs
^^^^^^

In Windows 7, ``conda install -c conda-forge nodejs`` might throw an
``IOError``. The workaround for this is to install ``yarn`` and ``nodejs``
outside of ``conda`` from https://yarnpkg.com (the stable version) and
https://nodejs.org (the LTS version), respectively.
