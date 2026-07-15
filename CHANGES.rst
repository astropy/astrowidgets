1.0 (unreleased)
================

TBD -- First stable release.

Bug Fixes
---------

- Fixed the bqplot ``ImageWidget`` so that mouse clicks no longer raise
  ``AttributeError`` and no longer block user-registered ``on_msg``
  callbacks. The dead click-to-center and interactive-marking code left
  over from before the switch to ``astro-image-display-api`` was removed;
  those features will be rebuilt (see #201). [#206]
