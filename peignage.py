#!/usr/bin/python
import math
import random
from gimpfu import *


def new_fiber_channels(img, fparams):
    """
    Return an array of correlated fiber channels.
    The length of the array equals the number of layers in img.
    """
    # Define constants
    w = img.width
    h = img.height
    l = fparams["length"]
    d = fparams["width"]
    r = fparams["width"] / 2.0
    # Internal frame
    w_in = img.width - 2*r
    h_in = img.height - 2*r
    # Trigonometric shortcuts
    cos_t = math.cos(math.radians(fparams["orientation"]))
    sin_t = math.sin(math.radians(fparams["orientation"]))
    # Grid size
    w_grid = h_in*sin_t + w_in*cos_t + l
    h_grid = h_in*cos_t + w_in*sin_t
    # Grid translation
    x_shift = r - l*cos_t + w_in*sin_t*sin_t
    y_shift = r - l*sin_t - w_in*sin_t*cos_t

    # Create base grid
    # Insert Intelligent Stuff Here


    # Set background and foreground colors
    gimp.set_background("black")
    gimp.set_foreground("white")

    # Instantiate black layers
    flayers = []
    for i in range(len(img.layers)):
        flayer = gimp.Layer(img, "Fibers", img.width, img.height,
                            RGBA_IMAGE, 100, LAYER_MODE_NORMAL)
        img.add_layer(flayer, 0)
        pdb.gimp_drawable_edit_fill(flayer, FILL_BACKGROUND)
        flayers.append(flayer)

    # Set base brush parameters
    pdb.gimp_context_set_brush("2. Hardness 050")
    pdb.gimp_context_set_brush_size(fparams["width"])
    pdb.gimp_context_set_brush_hardness(fparams["hardness"] / 10.0)

    # Paint over multiple layers
    # Insert Intelligent Stuff Here
    #pdb.gimp_paintbrush_default(flayers[rand], 4, [Ax, Ay, Bx, By])
    #pdb.gimp_paintbrush_default(flayers[0], 4, [100, 100, 300, 500])


    # Create channels based on brightness levels
    # White noise means full effect, black noise means no effect
    fchannels = []
    for flayer in flayers:
        pdb.gimp_selection_all(img)
        pdb.plug_in_colortoalpha(img, flayer, "black")
        pdb.gimp_image_select_item(img, CHANNEL_OP_REPLACE, flayer)
        img.remove_layer(flayer)
        fchannel = pdb.gimp_selection_save(img)
        fchannels.append(fchannel)

    return fchannels


def peignage(img,
             drawable,
             fibers_width=15,
             fibers_hardness=8,
             fibers_density=5,
             fibers_length=5,
             fibers_orientation=0,
             fibers_variation=0,
             overlap_gain=0):

    fparams = { "width": fibers_width,           # brush size
                "hardness": fibers_hardness,     # brush hardness
                "density": fibers_density,
                "length": fibers_length,
                "orientation": fibers_orientation,
                "variation": fibers_variation }

    #pdb.gimp_image_undo_group_start(img)

    # Store current background color
    bg_color_tmp = gimp.get_background()
    fg_color_tmp = gimp.get_foreground()

    # Create as many fiber channels as there are layers
    fiber_channels = new_fiber_channels(img, fparams)

    # Apply one fiber mask on each layer
    for i, fiber_channel in enumerate(fiber_channels):
        layer = img.layers[i]
        fiber_channel.name = layer.name
        pdb.gimp_image_set_active_channel(img, fiber_channel)
        fiber_mask = pdb.gimp_layer_create_mask(layer, ADD_MASK_CHANNEL)
        pdb.gimp_layer_add_mask(layer, fiber_mask)

    # Add brightness on intersecting regions
    if overlap_gain > 0:
        for layer in img.layers()[1:]:
            overlap = layer.copy()
            overlap.mode = ADDITION_MODE
            pdb.gimp_image_raise_item_to_top(img, overlap)
            pdb.gimp_drawable_brightness_contrast(overlap, overlap_gain / 20.0, 0)
            pdb.gimp_image_merge_down(img, overlap, CLIP_TO_IMAGE)

    # Merge all fiber layers
    pdb.gimp_image_merge_visible_layers(img, CLIP_TO_IMAGE)

    # Clean up
    pdb.gimp_selection_none(img)
    gimp.set_background(bg_color_tmp)
    gimp.set_foreground(fg_color_tmp)

    #pdb.gimp_image_undo_group_end(img)



    """
    # Extract islands
    pdb.gimp_context_set_antialias(False)
    pdb.gimp_context_set_feather(noise_smoothness)
    if noise_smoothness:
        pdb.gimp_context_set_feather_radius(50, 50)
    pdb.gimp_context_set_sample_threshold(0.5)
    pdb.gimp_image_select_color(img, CHANNEL_OP_REPLACE, layer_noise, "white")

    pdb.gimp_selection_invert(img)
    """


register(
 "python-fu-peignage",
 "Mingle layers like silk fibers",
 "Mingle layers like silk fibers",
 "Oriane Tury",
 "GNU GPLv3",
 "2021",
 "<Image>/Filters/Peignage",
 "RGB*",
 [
  (PF_SLIDER, "fibers_width", "Fibers Width", 15, (1, 1000, 1)),
  (PF_SLIDER, "fibers_hardness", "Fibers Hardness", 8, (0, 10, 1)),
  (PF_SLIDER, "fibers_density", "Fibers Density", 5, (1, 10, 1)),
  (PF_SLIDER, "fibers_length", "Fibers Length", 5, (1, 10, 1)),
  (PF_SLIDER, "fibers_orientation", "Fibers Orientation", 0, (-90, 90, 1)),
  (PF_SLIDER, "fibers_variation", "Fibers Variation", 0, (0, 10, 1)),
  (PF_SLIDER, "overlap_gain", "Overlap Gain", 0, (0, 10, 1)),
 ],
 [],
 peignage)

main()
