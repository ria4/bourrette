#!/usr/bin/python
import math
import random
from gimpfu import *

# Local import
import grid


# Hey! It's a debug switch!
DEBUG = False

# Constants for exponential input conversion
L_A = -0.000344125449991
L_B = 0.11153346297
L_C = -7.0189446165
W_A = 3.27285167629e-05
W_B = 0.0594682389829
W_C = -6.96725624648

# Mapping for density conversion
# Note that density_map[9] may take a *long* time
density_map = [.1, .25, .4, .5, .6, .7, .8, .9, .95, .99, .996]

# Constants for fiber drawing
LEFTOVERS_THRESHOLD = 100
COVERAGE_FAILSAFE = 100
FIBER_STROKING = False
GOTTA_GO_FAST = True
TRANSPARENCY_CHECK = False

# If the plug-in does not respond, you might want to use:
# sudo kill $(ps aux | grep 'peignage.py -gimp' | awk '{print $2}')


def get_coverage(img, drawable, max_area):
    # Return the ratio of (alpha-weighted) pixels over the padded frame
    pdb.gimp_selection_all(img)
    mean, std_dev, median, pixels, count, percentile = \
            pdb.gimp_drawable_histogram(drawable, HISTOGRAM_VALUE, 0, 1)
    return pixels * 1.0 / max_area


def apply_mask_and_merge(img, layer_base, layer_mask, layer_collage, soft_overlap):

    # Create channel based on brightness levels
    # White noise means full effect, black noise means no effect
    # We cannot use gimp_image_select_color because it does not preserve smoothness
    pdb.plug_in_colortoalpha(img, layer_mask, "black")
    pdb.gimp_image_select_item(img, CHANNEL_OP_REPLACE, layer_mask)
    channel_mask = pdb.gimp_selection_save(img)

    # Create buffer layer
    # It needs to be visible to be merged down, because reasons
    layer_buffer = layer_base.copy()
    layer_buffer.name = "Fibers Buffer"
    pdb.gimp_item_set_visible(layer_buffer, True)
    img.add_layer(layer_buffer, 1)

    # Add mask channel to buffer layer
    pdb.gimp_image_set_active_channel(img, channel_mask)
    mask_buffer = pdb.gimp_layer_create_mask(layer_buffer, ADD_MASK_CHANNEL)
    pdb.gimp_layer_add_mask(layer_buffer, mask_buffer)

    # Blend intersecting areas with softlight
    if soft_overlap:
        layer_buffer_overlap = layer_buffer.copy()
        layer_buffer_overlap.mode = LAYER_MODE_SOFTLIGHT
        img.add_layer(layer_buffer_overlap, 0)
        pdb.gimp_image_merge_down(img, layer_buffer_overlap, CLIP_TO_IMAGE)
        layer_collage = img.layers[0]

    # Merge collage layer into buffer layer
    pdb.gimp_image_merge_down(img, layer_collage, CLIP_TO_IMAGE)
    img.layers[0].name = "Fibers Collage"

    # Stroke lines around fibers (for demo purposes only)
    if FIBER_STROKING:
        gimp.set_foreground("black")
        pdb.gimp_drawable_edit_stroke_selection(img.layers[0])
        gimp.set_foreground("white")

    # Refresh display to show progress
    pdb.gimp_displays_flush()
    img.remove_channel(channel_mask)


def make_fiber_collage(img, fparams):

    # Insulate base layers
    n_layers = len(img.layers)
    base_layers = list(img.layers)
    for layer in base_layers:
        pdb.gimp_item_set_visible(layer, False)

    # Set background and foreground colors
    gimp.set_background("black")
    gimp.set_foreground("white")

    # Set base brush parameters
    pdb.gimp_context_set_brush("2. Hardness 050")
    pdb.gimp_context_set_brush_size(fparams["width"])
    pdb.gimp_context_set_brush_hardness(fparams["hardness"] / 10.0)

    # Set stroke line parameters
    if FIBER_STROKING:
        pdb.gimp_context_set_line_width(1)
        pdb.gimp_context_set_stroke_method(STROKE_LINE)

    # Instantiate black mask layer
    layer_mask = gimp.Layer(img, "Fibers Mask", img.width, img.height,
                            RGBA_IMAGE, 100, LAYER_MODE_NORMAL)
    pdb.gimp_item_set_visible(layer_mask, False)
    pdb.gimp_displays_flush()
    img.add_layer(layer_mask, n_layers)

    # Instantiate transparent collage layer
    layer_collage = gimp.Layer(img, "Fibers Collage", img.width, img.height,
                               RGBA_IMAGE, 100, LAYER_MODE_NORMAL)
    pdb.gimp_drawable_fill(layer_collage, FILL_TRANSPARENT)
    img.add_layer(layer_collage, 0)

    # Initialize grid
    g = grid.Grid(img, fparams)

    # Loop until non-transparent area is large enough
    coverage_loop_cnt = 0
    coverage = 0
    coverage_target = density_map[int(round(fparams["density"]))]

    # Very very rough prediction of first-round grid density
    decimate_factor = 0
    if coverage_target < .7:
        coverage_predict = 1 + math.pi * fparams["width"] / 4 / fparams["length"]
        coverage_predict *= fparams["hardness"] / 10.0 / grid.FIBER_INLINE_GAP
        coverage_predict = max(0, coverage_predict - .05)
        decimate_factor = 3.0/5 * coverage_target / coverage_predict

    while coverage < coverage_target and coverage_loop_cnt < COVERAGE_FAILSAFE:

        # Fill grid with proper fibers
        g.setup_fibers()

        # Lower density if necessary
        if decimate_factor:
            g.decimate(decimate_factor)

        # Pick speed
        if GOTTA_GO_FAST:

            # Compute number of fibers per buffer
            n_buffer = len(g.fibers) / n_layers

            # Pick a random base layer
            base_layer = random.choice(base_layers)
            random.shuffle(base_layers)

            # Loop to add fibers layer after layer
            for i in range(n_layers):

                # Register current base layer
                layer_base = base_layers[i]
                layer_collage = img.layers[0]

                # Reset mask layer to full black
                pdb.gimp_selection_all(img)
                pdb.gimp_drawable_fill(layer_mask, FILL_BACKGROUND)

                # Draw fibers on layer mask
                for j in range(n_buffer):
                    fiber = g.fibers[i*n_buffer + j]
                    coords = [fiber.m.x, fiber.m.y, fiber.n.x, fiber.n.y]
                    if TRANSPARENCY_CHECK:
                        # Remove partially transparent fibers (e.g. after reframing)
                        if fiber.has_transparent_end(layer_base):
                            continue
                    pdb.gimp_paintbrush_default(layer_mask, 4, coords)

                # Apply layer mask and merge to collage layer
                apply_mask_and_merge(img, layer_base, layer_mask,
                                     layer_collage, fparams["overlap"])

                # Update coverage to continue or break loop
                # Method _progress_set_text does not work so we use _init
                coverage = get_coverage(img, img.layers[0], g.max_area)
                pdb.gimp_progress_init('{0:.3g} / {1}'.format(coverage, coverage_target), None)
                if coverage > coverage_target:
                    break

        # Pick accuracy (in balancing between layers)
        else:

            # Compute number of fibers per buffer
            n_buffer = len(g.fibers) / n_layers**2
            buffers_end = n_buffer * n_layers**2
            n_leftovers = len(g.fibers) - buffers_end

            # Sometimes there may be only leftovers
            if n_buffer > 0:

                # Loop to prevent layer prioritization
                for i in range(n_layers):

                    # Cycle base layers
                    base_layers = base_layers[1:] + base_layers[:1]

                    # Loop to add fibers layer after layer
                    for j in range(n_layers):

                        # Register current base layer
                        layer_base = base_layers[j]
                        layer_collage = img.layers[0]

                        # Reset mask layer to full black
                        pdb.gimp_selection_all(img)
                        pdb.gimp_drawable_fill(layer_mask, FILL_BACKGROUND)

                        # Draw fibers on layer mask
                        for k in range(n_buffer):
                            fiber = g.fibers[(i*n_layers + j)*n_buffer + k]
                            coords = [fiber.m.x, fiber.m.y, fiber.n.x, fiber.n.y]
                            if TRANSPARENCY_CHECK:
                                if fiber.has_transparent_end(layer_base):
                                    continue
                            pdb.gimp_paintbrush_default(layer_mask, 4, coords)

                        # Apply layer mask and merge to collage layer
                        apply_mask_and_merge(img, layer_base, layer_mask,
                                             layer_collage, fparams["overlap"])

                        # Update coverage to continue or break loop
                        coverage = get_coverage(img, img.layers[0], g.max_area)
                        pdb.gimp_progress_init('{0:.3g} / {1}'.format(coverage, coverage_target), None)
                        if coverage > coverage_target:
                            break
                    if coverage > coverage_target:
                        break

            # Add leftover fibers when it makes a difference
            if len(g.fibers) < LEFTOVERS_THRESHOLD:

                # Loop to prevent prioritization
                for i in range(n_leftovers):

                    # Register current base layer
                    layer_base = base_layers[i % n_layers]
                    layer_collage = img.layers[0]

                    # Reset mask layer to full black
                    pdb.gimp_selection_all(img)
                    pdb.gimp_drawable_fill(layer_mask, FILL_BACKGROUND)

                    # Draw fiber on layer mask
                    fiber = g.fibers[buffers_end + i]
                    coords = [fiber.m.x, fiber.m.y, fiber.n.x, fiber.n.y]
                    if TRANSPARENCY_CHECK:
                        if fiber.has_transparent_end(layer_base):
                            continue
                    pdb.gimp_paintbrush_default(layer_mask, 4, coords)

                    # Apply layer mask and merge to collage layer
                    apply_mask_and_merge(img, layer_base, layer_mask,
                                         layer_collage, fparams["overlap"])

                    # Update coverage to continue or break loop
                    coverage = get_coverage(img, img.layers[0], g.max_area)
                    pdb.gimp_progress_init('{0:.3g} / {1}'.format(coverage, coverage_target), None)
                    if coverage > coverage_target:
                        break

        # Increment failsafe counter
        coverage_loop_cnt += 1

    # Clean up
    img.remove_layer(layer_mask)


def peignage(img_active,
             drawable_active,
             f_length=50,
             f_width=50,
             f_orientation=0,
             f_messiness=0,
             f_density=6,
             f_hardness=8,
             f_overlap=0):

    # Scale user input exponentially
    f_width_exp = img_active.height * math.exp(W_A*f_width**2 + W_B*f_width + W_C)
    f_length_exp = img_active.width * math.exp(L_A*f_length**2 + L_B*f_length + L_C)

    # Store all fiber parameters
    fparams = { "length": f_length_exp,
                "width": f_width_exp,
                "orientation": f_orientation,
                "messiness": f_messiness,
                "density": f_density,
                "hardness": f_hardness,
                "overlap": f_overlap }

    # Switch to a duplicate for destructive work
    img = pdb.gimp_image_duplicate(img_active)
    display = pdb.gimp_display_new(img)

    # Prevent overlarge undo stack
    if not DEBUG:
        pdb.gimp_image_undo_disable(img)
        #pdb.gimp_image_undo_group_start(img)

    # Store current background color
    bg_color_tmp = gimp.get_background()
    fg_color_tmp = gimp.get_foreground()

    # Extract fibers from all image layers
    make_fiber_collage(img, fparams)

    # Remove image layers
    for layer in img.layers[1:]:
        img.remove_layer(layer)

    # Clean up
    pdb.gimp_selection_none(img)
    gimp.set_background(bg_color_tmp)
    gimp.set_foreground(fg_color_tmp)

    if not DEBUG:
        pdb.gimp_image_undo_enable(img)
        #pdb.gimp_image_undo_group_end(img)


register(
 "python-fu-peignage",
 "Mingle layers like silk waste fibers",
 "Mingle layers like silk waste fibers",
 "Oriane Tury",
 "GNU GPLv3",
 "2021",
 "<Image>/Filters/Artistic/Peignage",
 "RGB*",
 [
  (PF_SLIDER, "f_length", "Fibers Length", 50, (1, 100, 1)),
  (PF_SLIDER, "f_width", "Fibers Width", 50, (1, 100, 1)),
  (PF_SLIDER, "f_orientation", "Fibers Orientation", 0, (-90, 90, 1)),
  (PF_SLIDER, "f_messiness", "Fibers Messiness", 0, (0, 10, 1)),
  (PF_SLIDER, "f_density", "Fibers Density", 6, (0, 10, 1)),
  (PF_SLIDER, "f_hardness", "Fibers Hardness", 8, (0, 10, 1)),
  (PF_SLIDER, "f_overlap", "Soft Overlap", 0, (0, 1, 1)),
  # i don't care much for the bulky PF_BOOL button
 ],
 [],
 peignage)

main()
