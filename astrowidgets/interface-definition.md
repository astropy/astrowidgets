# Image widget interface definition
## Intended to make discussion of API separate from implementation

## An astro image viewer has these attributes and methods and follows these conventions

### Conventions

1. The pixel at index `[0, 0]` is display in the lower left by default.
1. The center of a pixel is 0, i.e. the first pixel extends from `-0.5` to `0.5` -- the thing that needs to be addressed here is that APE-14 WCS makes statements about what it uses as pixel center and we should address that here.
1. ...?

### Attributes

Unless otherwise noted, all are settable.

1. `viewer` -- object that is doing the actual display on the screen. **READ ONLY**
1. `image_height` -- height of image on screen, in pixels
1. `image_width` -- width of image on screen, in pixels
1. ~~`pixel_offset`~~ -- offset added to each pixel before displaying it. **READ ONLY** **DISCUSSION item:** What is the display convention? Is the center of the "first" pixel at location (0, 0) or is it (0.5, 0.5)? In other words, is the center of a pixel 0 or 0.5?
1. `zoom_level` -- current zoom factor of the image, with the convention that `zoom_level=1` means one pixel in the image is 1 pixel on the screen.
1. `is_marking` -- `True` if interactive marking is in process, `False` otherwise. Setting to `True` has the side effect of disabling some other interactions.
1. `marker_style` -- the style currently being used to make marks. *Note:* Name change from `marker` to `marker_style` proposed in #145.
1. `stretch_options` -- a list of the options for stretching the data. **DISCUSSION topic:** Is this a list of strings? Astropy stretch objects? Something else?
1. `stretch` -- current image stretch. **DISCUSSION topic:** What are the values allowed here? String (e.g. `sqrt`) or astropy stretch object or either?
1. `autocut_options` -- list of options for auto-cutting the data. **DISCUSSION topic:** Is this a list of strings? Astropy stretch objects? Something else? **DISCUSSION topic:** Should these be called intervals or should the documentation of it at least mention that they are called intervals in the astropy documentation?
1. `cuts` -- current cuts. **DISCUSSION topic:** What are the values allowed here? `tuple` (e.g. `(2, 1000)`) or astropy interval object or either?
1. `cursor` -- location of cursor display. One of `top` or `bottom` or `None`.
1. `click_center` -- if `True`, clicking on image re-centers image on that location.
1. `click_drag` -- if `True`, pan the image by clicking and dragging.
1. `scroll_pan` -- if `True`, pan the image by scrolling.

### Methods

1. `load_fits(self, filename, **kwargs)` -- load a fits file into the viewer.
1. `load_nddata(self, nddata, **kwargs)` -- load an nddata object into viewer.
1. `def load_array(self, arr, **kwargs)` -- load an array into viewer.
1. `center_on(self, point)` -- center the viewer on `point` which may be a pixel or a `SkyCoord`.
1. `offset_by(self, dx, dy)` -- amount by which to offset the viewer. After applying the offset, the new center is (old center 1. offset). `dx` and `dy` can be either pixels or an angle quantity. **DISCUSSION item:** Why separate dx, dy instead of a point?
1. `zoom(self, value)` -- factor by which to zoom in or out, i.e. the factor by which to change `zoom_level` from its current value. **DISCUSSION item:** Would `zoom_by` or `change_zoom_by` be a better name?
1. `start_marking(self, marker_name=None, marker_style=None)` -- activate marking mode with option name and/or style for the markers. *Note:* Keyword argument name change to `marker_style` proposed in #145.
1. `stop_marking(self, clear_markers=False)` -- deactivate marking mode, optionally clearing the markers away.
1. `get_marker_names(self)` -- get list of the marker names currently in use.
1. `get_markers_by_name(self, marker_name, x_colname='x', y_colname='y', skycoord_colname='coord')` -- get information about markers with name `marker_name`. **DISCUSSION item:** The values for x, y depend on the convention adopted for pixel center.
1. `get_all_markers(self, x_colname='x', y_colname='y', skycoord_colname='coord')` -- get all of the markers. **DISCUSSION item:** Does the return table include marker names?
1. `add_markers(self, table, x_colname='x', y_colname='y', skycoord_colname='coord', use_skycoord=False, marker_name=None)` -- add markers from an astropy table programmatically. **DISCUSSION itme:** Shouldn't this have a `marker_style` option too?
1. `remove_markers_by_name(self, marker_name)` -- remove markers with name `marker_name`
1. `remove_all_markers(self)` -- remove all of the markers.
1. `validate_marker_name(self, marker_name)` -- checks whether marker name is used internally.
1. `set_cached_state(self)` -- cache several settings before starting interactive marking.
1. `restore_and_clear_cached_state(self)` -- restore settings when interactive marking ends.
1. `save(self, filename, overwrite=False)` -- save image, with format determined by the extension of the `filename`. **DISCUSSION item:** Hmmm, wouldn't `clobber` be a more intuitive name than `overwrite`? jk, jk.....
