#!/usr/bin/python
import random
from gimpfu import *


def new_fiber_channel(img, fparams):
    # Set background and foreground colors
    gimp.set_background("black")
    gimp.set_foreground("white")

    # Instantiate black layer
    flayer = gimp.Layer(img, "Fibers", img.width, img.height,
                        RGBA_IMAGE, 100, LAYER_MODE_NORMAL)
    img.add_layer(flayer, 0)
    pdb.gimp_drawable_edit_fill(flayer, FILL_BACKGROUND)

    # Create channel based on brightness levels
    # White noise means full effect, black noise means no effect
    pdb.gimp_selection_all(img)
    pdb.plug_in_colortoalpha(img, flayer, "black")
    pdb.gimp_image_select_item(img, CHANNEL_OP_REPLACE, flayer)
    img.remove_layer(flayer)
    fchannel = pdb.gimp_selection_save(img)

    return fchannel


def peignage(img,
             drawable,
             fibers_width=15,
             fibers_hardness=5,
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

    # Apply a different fiber mask on every layer
    for layer in img.layers:
        fiber_channel = new_fiber_channel(img, fparams)
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
  (PF_SLIDER, "fibers_hardness", "Fibers Hardness", 5, (0, 10, 1)),
  (PF_SLIDER, "fibers_density", "Fibers Density", 5, (1, 10, 1)),
  (PF_SLIDER, "fibers_length", "Fibers Length", 5, (1, 10, 1)),
  (PF_SLIDER, "fibers_orientation", "Fibers Orientation", 0, (-90, 90, 1)),
  (PF_SLIDER, "fibers_variation", "Fibers Variation", 0, (0, 10, 1)),
  (PF_SLIDER, "overlap_gain", "Overlap Gain", 0, (0, 10, 1)),
 ],
 [],
 peignage)

main()
