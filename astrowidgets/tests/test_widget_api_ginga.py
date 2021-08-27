import pytest

from .widget_api_test import ImageWidgetAPITest

ginga = pytest.importorskip("ginga",
                            reason="Package required for test is not "
                                   "available.")
from astrowidgets.ginga import ImageWidget  # noqa: E402


class TestGingaWidget(ImageWidgetAPITest):
       def setup_class(self):
              self.image_widget_class = ImageWidget
              super().setup_class(self)
