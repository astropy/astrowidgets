1.0 (unreleased)
================

TBD -- First stable release.

New Features
------------

- The cursor-information readout (X/Y, RA/Dec, pixel value) is now shared
  between the ginga and bqplot backends via
  ``astrowidgets.cursor_info.CursorInfoMixin``. The ``cursor`` property
  ('top', 'bottom', or ``None``), which had been dropped from the bqplot
  backend, is available again on both backends, and a new
  ``sky_coordinate_format`` property selects between decimal degrees (the
  new default for both backends) and sexagesimal for the RA/Dec display.
  As part of this change the ginga readout switched from sexagesimal to
  decimal degrees by default; set ``sky_coordinate_format`` to
  ``'sexagesimal'`` for the previous behavior.

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
