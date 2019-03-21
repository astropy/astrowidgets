"""Module containing some functionalities of ``astrowidgets``."""
from __future__ import print_function, unicode_literals

# STDLIB
import sys
import os
import random
import functools
import warnings
import logging
from traitlets import Unicode, validate

# THIRD-PARTY
import numpy as np
import astropy
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack
import astropy.visualization as aviz
from astropy.nddata import NDData
from astropy.nddata.utils import block_reduce
from astropy.wcs import WCS

# Jupyter widgets
import ipywidgets as widgets
from ipywidgets import interact, interactive, fixed, interact_manual

# IPython
from IPython.display import display

# Matplotlib
import matplotlib
# Make sure that we are usign TkAgg
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import matplotlib.widgets as matwidgets 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.backend_tools import SetCursorBase, ToolPan

__all__ = ['MatImageWidget']


# Allowed locations for cursor display
ALLOWED_CURSOR_LOCATIONS = ['top', 'bottom', None]

# List of marker names that are for internal use only
RESERVED_MARKER_SET_NAMES = ['all']

import tkinter as Tk

class MatImageWidget(widgets.DOMWidget):
    """
    Image widget for Jupyter notebook using Matplotlib as a non-interactive backend.

    .. todo:: Any property passed to constructor has to be valid keyword.

    Parameters
    ----------
    Please fill in the blanks...

    """
    #constructor
    def __init__(self, width =8, height = 6, pixel_coords_offstet=0, *args, **kwargs):
        super().__init__()
        
        _view_name = Unicode('astrowidgets').tag(sync=True)
        _view_module = Unicode('Matplotlib').tag(sync=True)
        _view_module_version = Unicode('0.1.0').tag(sync=True)

        self._jup_img = widgets.Image(format='png')

        # Set the image margin to over the widgets default of 2px on all sides.
        self._jup_img.layout.margin = '0'

        # Set both of those to ensure consistent display in notebook
        # and jupyterlab when the image is put into a container smaller
        # than the image.

        self._jup_img.max_width = '100%'
        self._jup_img.height = '100%'

        # Set the width of the box containing the image to the desired width
        #self.layout.width = str(image_width)

        # Marker
        self.marker = {'alpha':'0.25', 'marker': 'o', 'color': 'blue', 's': 20}

        # Maintain marker tags as a set because we do not want
        # duplicate names.
        self._marktags = set()
        # Let's have a default name for the tag too:
        self._default_mark_tag_name = 'default-marker-name'
        self._interactive_marker_set_name = 'interactive-markers'

        # coordinates display
        self._jup_coord = widgets.HTML('Coordinates show up here')
        # This needs ipyevents 0.3.1 to work
        #self.draw_cursor(event) = canvas.draw_cursor(self, event)

        # Define a callback that shows the output of a print
        self._print_out = widgets.Output()

        self._cursor = 'bottom'
        self.children = [self._jup_img, self._jup_coord]


        # set DEBUG for everything
        logging.basicConfig(level=logging.DEBUG)
        logger = logging.getLogger('matplotlib')

        # set WARNING for Matplotlib
        logger.setLevel(logging.WARNING)


#    def _repr_html_(self):
#        """
#        Show widget in Jupyter notebook.
#        """
#        from IPython.display import display
#        return display(self._widget)

    def load_fits(self, fitsorfn, numhdu, memmap=None):
        """
        Loads a FITS file into the viewer.

        Parameters
        ----------
        fitsorfn : str or HDU
            Either a file name or an HDU (*not* an HDUList).
            If file name is given, WCS in primary header is automatically
            inherited. If a single HDU is given, WCS must be in the HDU
            header.

        numhdu : int or ``None``
            Extension number of the desired HDU.
            If ``None``, it is determined automatically.

        memmap : bool or ``None``
            Memory mapping.
            If ``None``, it is determined automatically.

        """
        # Will need to tweak this a bit to make it work.
        if isinstance(fitsorfn, str):
            hdul = fits.open(fitsorfn)
            image = hdul[numhdu].data
            hdr = hdul[numhdu].header
        elif isinstance(fitsorfn, (fits.ImageHDU, fits.CompImageHDU, fits.PrimaryHDU)):
            hdu = fits.open(fitsorfn)
            image = hdu[0].data
            hdr = hdu[0].header
        return image, hdr

    def load_nddata(self, ndd):  # Wondering how come I cannot put a ``self`` in here before ``ndd``? 
        """
        Load an ``NDData`` object into the viewer.

        Parameters
        ----------
        ndd : `~astropy.nddata.NDData`
            ``NDData`` with image data and WCS.

        """
        if isinstance(ndd, NDData):
            image = ndd.data
            if ndd.wcs:    
                wcs = ndd.wcs
                hdr = wcs.to_header()
                return image, hdr
            else: 
                raise ValueError("No basic FITS header info provided.")
        
    #def load_array(self, arr):
        """
        Load a 2D array into the viewer.
        .. note:: Use :meth:`load_nddata` for WCS support.
        Parameters
        ----------
        arr : array-like
            2D array.
        """
        #if isinstance(arr, (list, tuple)):
        #    image = np.array(arr)
        #elif isinstance(arr, np.ndarray):
        #    image = arr
        #hdr = None
        #return image


#    # Will need to tweak this a bit to make it work in this context.
#    def plot_click(self, image):
#        """Plot dummy data and handle user clicks."""
#
#        # Show data
#        if image is not False:
#            plt.imshow(image)
#
#        # Will need to tweak this a bit to make it work in this context.
#        def _on_click(event):
#            """Print and mark clicked coordinates."""
#
#            # Print data coordinates to screen
#            print('X, Y:', event.xdata, event.ydata)
#
#            # Mark them on plot too
#            plt.plot(event.xdata, event.ydata, 'rx')
#
#        # Register click events
#        plt.connect('button_press_event', _on_click)


    def show_image(self, image, hdr, percl=99, percu=None, figsize=(8, 6),
               cmap='viridis', log=True,
               show_colorbar=True, show_ticks=True,
               fig=None, ax=None, input_ratio=None):
        """
        Show an image in matplotlib with some basic astronomically-appropriate stretching.
  
        """
        #root = Tk.Tk()
        #root.wm_title("Embedding in TK")

        wcs = WCS(hdr)

        if percu is None:
            percu = percl
            percl = 100 - percl

        if (fig is None and ax is not None) or (fig is not None and ax is None):
            raise ValueError('Must provide both "fig" and "ax" '
                             'if you provide one of them')
        elif fig is None and ax is None:
            # fig, ax = plt.subplots(1, 1, figsize=figsize)
            fig = plt.figure(figsize=figsize)
            ax = plt.subplot(projection=wcs)
            if figsize is not None:
                # Rescale the fig size to match the image dimensions, roughly
                image_aspect_ratio = image.shape[0] / image.shape[1]
                figsize = (max(figsize) * image_aspect_ratio, max(figsize))
                #print(figsize)

        # To preserve details we should *really* downsample correctly and
        # not rely on matplotlib to do it correctly for us (it won't).

        # So, calculate the size of the figure in pixels, block_reduce to
        # roughly that, and display the block-reduced image.
        fig_size_pix = fig.get_size_inches() * fig.dpi

        ratio = (image.shape // fig_size_pix).max()

        if ratio < 1:
            ratio = 1

        ratio = input_ratio or ratio

        # Divide by the square of the ratio to keep the flux the same in the
        # reduced image
        reduced_data = block_reduce(image, ratio) / ratio**2

        # Of course, now that we have down-sampled, the axis limits are changed to
        # match the smaller image size. Setting the extent will do the trick to
        # change the axis display back to showing the actual extent of the image.
        extent = [0, image.shape[1], 0, image.shape[0]]

        if log:
            stretch = aviz.LogStretch()
        else:
            stretch = aviz.LinearStretch()

        norm = aviz.ImageNormalize(reduced_data, interval=aviz.AsymmetricPercentileInterval(percl, percu),
                                      stretch=stretch)
        
        # The following line makes it so that the zoom level no longer changes,
        # otherwise Matplotlib has a tendency to zoom out when adding overlays.
        # ax.set_autoscale_on(False)

        # Can also try to make the colorbar the same height as the image by setting
        # aspect='auto'. 
        im = ax.imshow(reduced_data, norm=norm, origin='lower',
                       cmap=cmap, extent=extent, aspect='equal')
        
        plt.xlabel(r'$RA$')
        plt.ylabel(r'$Dec$')

        if show_colorbar:
            # I haven't a clue why the fraction and pad arguments below work to make
            # the colorbar the same height as the image, but they do....unless the image
            # is wider than it is tall. Sticking with this for now anyway...
            # Thanks: https://stackoverflow.com/a/26720422/3486425
            #fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

            # I have eliminated the need for using the ``fraction`` and ``pad`` attributes of 
            # the colorbar by simply adding the ``aspect`` attribute to imshow() instead. 
            fig.colorbar(im, ax=ax)

        # The plot now shows ticks at the bottom and left sides, while the colorbar occupies a
        # position on the right-hand side of the diagram. 
        if not show_ticks:
            ax.tick_params(labelbottom=False, labelleft=False)


        # GUI
        root = Tk.Tk()
        root.wm_title("Embedding in Tk")

        # a tk.DrawingAreaanvas = FigureCanvasT
        canvas = FigureCanvasTkAgg(fig, master=root)  # A Tk.DrawingArea
        canvas.draw()
        #canvas.get_tk_widget().grid(row=0, column=1)
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        #_canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
