#!/usr/bin/python
import random
from gimpfu import *


def decreusage(img,
               drawable,
               noise_seed=0,
               noise_detail=1,
               noise_size=9,
               noise_smoothness=False,
               noise_threshold=15,
               hue_shift=0,
               chroma_shift=30):

    pdb.gimp_image_undo_group_start(img)

    # Shift inputs to internal values
    noise_size = 11 - noise_size

    # Store current background color
    bg_color_tmp = gimp.get_background()

    # Register layer positon
    drawable_position = pdb.gimp_image_get_item_position(img, drawable)

    # Create noise layer
    layer_noise = gimp.Layer(img, "Noise %s" % drawable.name,
                             img.width, img.height,
                             RGBA_IMAGE, 100, LAYER_MODE_NORMAL)
    img.add_layer(layer_noise, drawable_position)
    if noise_seed == 0:
        noise_seed = random.getrandbits(32)
    pdb.plug_in_solid_noise(img, layer_noise, False, False,
                            noise_seed, noise_detail, noise_size, noise_size)

    # Extract islands
    pdb.gimp_context_set_antialias(False)
    pdb.gimp_context_set_feather(noise_smoothness)
    if noise_smoothness:
        pdb.gimp_context_set_feather_radius(50, 50)
    pdb.gimp_context_set_sample_threshold(0.5)
    pdb.gimp_image_select_color(img, CHANNEL_OP_REPLACE, layer_noise, "white")
    gimp.set_background("white")
    pdb.gimp_drawable_edit_fill(layer_noise, FILL_BACKGROUND)

    # Apply noise threshold
    threshold = int(noise_threshold * 255 / 100)
    pdb.gimp_selection_invert(img)
    gimp.set_background(threshold, threshold, threshold)
    pdb.gimp_drawable_edit_fill(layer_noise, FILL_BACKGROUND)

    # Create channel based on noise brightness levels
    # White noise means full effect, black noise means no effect
    pdb.gimp_selection_all(img)
    pdb.plug_in_colortoalpha(img, layer_noise, "black")
    pdb.gimp_image_select_item(img, CHANNEL_OP_REPLACE, layer_noise)
    channel_noise = pdb.gimp_selection_save(img)
    channel_noise.name = layer_noise.name

    # Add channel mask to drawable duplicate
    layer_dup = drawable.copy()
    layer_dup.name = "Duplicate %s" % drawable.name
    pdb.gimp_image_set_active_channel(img, channel_noise)
    mask_noise = pdb.gimp_layer_create_mask(layer_dup, ADD_MASK_CHANNEL)
    pdb.gimp_layer_add_mask(layer_dup, mask_noise)

    # Apply hue-chroma mask
    img.add_layer(layer_dup, drawable_position + 1)
    gegl_parameters = "hue-chroma hue=%s chroma=%s" % (hue_shift, chroma_shift)
    pdb.gegl_gegl(layer_dup, gegl_parameters)

    # Clean up
    img.remove_layer(layer_noise)
    pdb.gimp_image_merge_down(img, layer_dup, CLIP_TO_IMAGE)
    pdb.gimp_selection_none(img)
    gimp.set_background(bg_color_tmp)

    pdb.gimp_image_undo_group_end(img)


register(
 "python-fu-decreusage",
 "Chromatize layer over a solid noise pattern",
 "Chromatize layer over a solid noise pattern",
 "Oriane Tury",
 "GNU GPLv3",
 "2021",
 "<Image>/Filters/Decreusage",
 "RGB*",
 [
  (PF_INT32, "noise_seed", "Noise Seed", 0),
  (PF_SPINNER, "noise_detail", "Noise Detail", 1, (0, 5, 1)),
  (PF_SLIDER, "noise_size", "Noise Size", 9, (1, 10, 1)),
  (PF_BOOL, "noise_smoothness", "Noise Smoothness", False),
  (PF_SLIDER, "noise_threshold", "Noise Threshold", 15, (0, 100, 1)),
  (PF_SLIDER, "hue_shift", "Hue Adjustment", 0, (-180, 180, 1)),
  (PF_SLIDER, "chroma_shift", "Chroma Adjustment", 30, (-180, 180, 1)),
 ],
 [],
 decreusage)

main()
