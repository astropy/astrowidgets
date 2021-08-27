import pytest

from .widget_api_test import ImageWidgetAPITest

_ = pytest.importorskip("bqplot",
                        reason="Package required for test is not "
                               "available.")
from astrowidgets.bqplot import ImageWidget  # noqa: E402


class TestGingaWidget(ImageWidgetAPITest):
       def setup_class(self):
              self.image_widget_class = ImageWidget
              super().setup_class(self)
