1.0 (unreleased)
================

TBD -- First stable release.

New Features
------------

- The cursor-information readout (X/Y, RA/Dec, pixel value) is now shared
  between the ginga and bqplot backends via
  ``astrowidgets.cursor_info.CursorInfoMixin``. Both backends have a
  ``cursor`` property (``'top'``, ``'bottom'``, or ``None``) and a new
  ``sky_coordinate_format`` property that selects decimal degrees (the
  default) or sexagesimal for the RA/Dec display. The ginga readout now
  defaults to decimal degrees; set ``sky_coordinate_format`` to
  ``'sexagesimal'`` for the previous display. [#222]

- The bqplot ``ImageWidget`` no longer sends the whole image to the
  browser each time the cuts, stretch or colormap change. The pixels are
  sent once, when an image is loaded, and the cuts and stretch are applied
  in the browser, through the limits and the colors of the color scale. A
  stretch too steep to be drawn that way without banding, such as the
  default ``LogStretch``, is still applied to the pixels and sent when the
  cuts or stretch change, but as 8 bit, a quarter of the size. Removing
  catalogs no longer sends the image either. [#244]

Bug Fixes
---------

- Fixed the RA/Dec shown in the bqplot ``ImageWidget`` cursor readout,
  which transposed the x and y pixel coordinates when converting to sky
  coordinates. Also fixed both backends reading pixel 0's value for
  cursor positions slightly outside the image on the negative side.

- Fixed the bqplot ``ImageWidget`` so that mouse clicks no longer raise
  ``AttributeError`` and no longer block user-registered ``on_msg``
  callbacks. The dead click-to-center and interactive-marking code left
  over from before the switch to ``astro-image-display-api`` was removed;
  those features will be rebuilt (see #201). [#206]
