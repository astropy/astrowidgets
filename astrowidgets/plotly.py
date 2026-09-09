import ipywidgets as ipw
from astro_image_display_api.image_viewer_logic import ImageViewerLogic
import astropy.visualization as apviz
import numpy as np
from astropy.nddata import NDData
import plotly.graph_objects as go
import plotly.express as px

def docs_from_super_if_missing(cls):
    """
    Decorator to copy the docstrings from the interface methods to the
    methods in the class.
    """
    for name, method in cls.__dict__.items():
        if not name.startswith("_"):
            if method.__doc__:
                continue
            interface_method = getattr(ImageViewerLogic, name, None)

            if interface_method:
                method.__doc__ = interface_method.__doc__
    return cls

@docs_from_super_if_missing
class ImageWidget(ipw.VBox, ImageViewerLogic):
    def __init__(self, *args, display_width=500, display_aspect_ratio=1):
        super().__init__(*args)
        print('Initializing ImageWidget')
        self._set_up_catalog_image_dicts()

        self._astro_im = go.FigureWidget(data=[go.Heatmap(z=np.zeros((10, 10)))])
        self._default_cuts = apviz.MinMaxInterval()
        self._default_stretch = None

        self._data = None
        self._wcs = None

        # Use this to manage whether or not to send changes in zoom level
        # to the viewer.
        self._viewport_change_source_is_gui = False

        # Provide an Output widget to which prints can be directed for
        # debugging.
        self._print_out = ipw.Output()

        self.marker = {'color': 'red', 'radius': 20, 'type': 'square'}
        self._cuts = apviz.AsymmetricPercentileInterval(1, 99)

        self._cursor = ipw.HTML('Coordinates show up here')

        self.children = [self._astro_im, self._cursor]

    def _interval_and_stretch(self, stretch=None, cuts=None):
        """
        Stretch and normalize the data before sending to the viewer.
        """
        interval = cuts or self._default_cuts
        intervaled = interval(np.asarray(self._data))

        stretch = stretch or self._default_stretch
        if stretch:
            stretched = stretch(intervaled)
        else:
            stretched = intervaled

        return stretched

    def _send_data(self, reset_view=True, stretch=None, cuts=None):
        if self._data is not None:
            print('Sending data to viewer')
            fig = px.imshow(
                self._interval_and_stretch(stretch=stretch, cuts=cuts)
            )
            print('Data sent to viewer', fig.data)
            self._astro_im.add_trace(fig.data[0])


    @property
    def _current_image_label(self):
        """
        Image label for the most recently loaded image
        """
        return list(self._images.keys())[-1]

    def set_stretch(self, value, image_label=None, **kwargs):
        super().set_stretch(value, image_label=image_label, **kwargs)
        self._send_data(stretch=value)

    def set_cuts(self, value, image_label=None, **kwargs):
        super().set_cuts(value, image_label=image_label, **kwargs)
        self._send_data(cuts=self.get_cuts(image_label=image_label))

    @property
    def viewer(self):
        return self._astro_im

    # The methods, grouped loosely by purpose
    def load_image(self, image, image_label=None, **kwargs):
        super().load_image(image, image_label=image_label, **kwargs)
        data = self.get_image(image_label=image_label)

        self._data = data.data if isinstance(data, NDData) else data
        self._send_data()

    def set_colormap(self, cmap_name, image_label=None, **kwargs):
        super().set_colormap(cmap_name, image_label=image_label, **kwargs)
        self._astro_im.set_color(bqcolors(cmap_name, reverse=False))

    def load_catalog(
        self,
        table,
        **kwargs
    ):
        super().load_catalog(table, **kwargs)
        catalog_label = kwargs.pop("catalog_label", None)
        this_catalog = self.get_catalog(catalog_label=catalog_label)
        self._astro_im.plot_named_markers(
            this_catalog["x"],
            this_catalog["y"],
            str(catalog_label),
            **self.get_catalog_style(catalog_label=catalog_label)
        )

    def set_catalog_style(
            self,
            catalog_label=None,
            shape="circle",
            color="red",
            size=5,
            **kwargs
        ):
        super().set_catalog_style(
            catalog_label=catalog_label,
            shape=shape,
            color=color,
            size=size,
            **kwargs
        )
        this_catalog = self.get_catalog(catalog_label=catalog_label)

        self._astro_im.plot_named_markers(
            this_catalog["x"],
            this_catalog["y"],
            str(catalog_label),
            color=color,
            size=size**2,  # bqplot expects size in pixels squared
            shape=shape,
        )

    @property
    def print_out(self):
        """
        Return an output widget for display in the notebook which
        captures any printed output produced by the viewer widget.

        Intended primarily for debugging.
        """
        return self._print_out
