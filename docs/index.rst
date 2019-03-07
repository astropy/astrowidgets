************
astrowidgets
************

``astrowidgets`` aims to be a set of astronomy widgets for Jupyter Lab
or Notebook, leveraging the Astropy ecosystem.

.. warning::

    ``astrowidgets`` is currently in very early development.
    Nothing is guaranteed to work or do anything in particular at this stage.


Installation
============

Currently, ``astrowidgets`` are only available via GitHub::

    pip install git+https://github.com/astropy/astrowidgets.git@master

.. note::

    If you are using a virtual environment and
    are experiencing problems with the installation, be sure to check and use
    the suitable pip version. Also, if you are having problem even with the following 
    dependencies installed, please make sure you have the most up-to-date versions. 
    Upgrade dependencies when and where necessary.

The following dependencies are needed *in the kernel of execution*
to use ``astrowidgets`` in either Jupyter Lab or Notebook:

* ``python >= 3.5``
* ``numpy``
* ``astropy``
* ``ipywidgets``
* ``ipyevents >= 0.3.1``

For Jupyter Notebook only:

* ``jupyter notebook``

For Jupyter Lab only:

* ``jupyter lab``
* ``nodejs``
* After installing dependencies, run::

    jupyter labextension install @jupyter-widgets/jupyterlab-manager ipyevents

For those using ``conda``, dependencies from the ``conda-forge`` channel
should be sufficient unless stated otherwise.


Widget with Ginga Toolkit
-------------------------

To use the widget with `Ginga <http://ginga.readthedocs.io>`_ toolkit,
you also need to install:

* ``ginga>=2.7.1``
* ``pillow``
* ``freetype``
* ``aggdraw``
* ``opencv`` (optional)

.. note::

    For vectorized drawing in ``aggdraw``, you can clone
    https://github.com/ejeschke/aggdraw/ and install its ``vectorized-drawing``
    branch from source.


Notes for Windows Users
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


Using astrowidgets
==================

.. toctree::
    :maxdepth: 2

    astrowidgets/index


Reference/API
=============

.. automodapi:: astrowidgets
    :no-main-docstr:
    :no-inheritance-diagram:
