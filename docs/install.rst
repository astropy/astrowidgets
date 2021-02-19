.. _astrowidgets_install:

Installation
============

This page contains the installation instructions for the abstract class in
``astrowidgets``. To use the concrete implementation with Ginga, please also
see :ref:`ginga_backend`.

.. _astrowidgets_install_pip:

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

.. _astrowidgets_install_conda:

Install with conda
------------------

``conda`` installation::

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

* ``python>=3.6``
* ``numpy``
* ``astropy``
* ``ipywidgets>=7.5``
* ``ipyevents>=0.6.3``
* ``jupyterlab>=1``
* ``nodejs``

After installing dependencies, for Jupyter Lab, run::

    jupyter labextension install @jupyter-widgets/jupyterlab-manager

For those using ``conda``, dependencies from the ``conda-forge`` channel
should be sufficient unless stated otherwise.

nodejs on Windows
^^^^^^^^^^^^^^^^^

In Windows 7, ``conda install -c conda-forge nodejs`` might throw an
``IOError``. The workaround for this is to install ``yarn`` and ``nodejs``
outside of ``conda`` from https://yarnpkg.com (the stable version) and
https://nodejs.org (the LTS version), respectively.
