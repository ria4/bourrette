#!/usr/bin/python
import math
import random
from gimpfu import *


# Constants for exponential input conversion
L_A = -0.000344125449991
L_B = 0.11153346297
L_C = -7.0189446165
W_A = 3.27285167629e-05
W_B = 0.0594682389829
W_C = -6.96725624648


class Grid:
    def __init__(self, img, fparams):
        # Imported constants
        self.w_img = img.width
        self.h_img = img.height
        self.l = fparams["length"]
        self.d = fparams["width"]
        self.r = fparams["width"] / 2.0
        # Internal frame
        self.w_in = self.w_img - 2*r
        self.h_in = self.h_img - 2*r
        # Precomputed trigonometrics
        self.cos_t = math.cos(math.radians(fparams["orientation"]))
        self.sin_t = math.sin(math.radians(fparams["orientation"]))
        # Dimensions
        self.w = self.h_in*self.sin_t + self.w_in*self.cos_t + self.l
        self.h = self.h_in*self.cos_t + self.w_in*self.sin_t
        # Cells
        self.n_lines = int(self.h / self.d) + 2
        self.grid_x_0 = [i * self.l / self.n_lines for i in range(self.n_lines)]
        self.n_columns = int(self.w / self.l)
        # Translation
        self.tx = self.r - self.l*self.cos_t + self.w_in*self.sin_t*self.sin_t
        self.ty = self.r - self.l*self.sin_t - self.w_in*self.sin_t*self.cos_t

        self.reset()

    def reset(self):
        self.coordinates = []
        grid_y_0 = -self.d * random.random()
        random.shuffle(self.grid_x_0)
        for i in range(self.n_lines):
            line_abscissae = []
            for j in range(self.n_columns):
                line_abscissae.append(self.grid_x_0[i] + j*self.l)
            self.coordinates.append([line_abscissae, grid_y_0 + i*self.d])

    def rotate_translate(self):
        coordinates = []
        for i in range(self.n_lines):
            y_old = self.coordinates[i][1]
            line_abscissae = []
            for j in range(self.n_columns):
                x_old = self.coordinates[i][0][j]
                x_new = x_old*self.cos_t - y_old*self.sin_t + self.tx
                grid_line_abscissae.append(x_new)
            y_new = x_old*self.sin_t + y_old*self.cos_t + self.ty
            coordinates.append([grid_line_abscissae, y_new])
        self.coordinates = coordinates

    # def creer points couples
    # tester
    # def projeter points aux frontières
    # def filtrer segments trop courts
    # def decimer


def new_fiber_channels(img, fparams):
    """
    Return an array of correlated fiber channels.
    The length of the array equals the number of layers in img.
    """
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
             f_length=50,
             f_width=50,
             f_orientation=0,
             f_variation=0,
             f_density=5,
             f_hardness=8,
             overlap_gain=0):

    # Scale user input exponentially
    f_width_exp = img.height * math.exp(W_A*f_width**2 + W_B*f_width + W_C)
    f_length_exp = img.width * math.exp(L_A*f_length**2 + L_B*f_length + L_C)

    # Store all fiber parameters
    fparams = { "length": f_length_exp,
                "width": f_width_exp,
                "orientation": f_orientation,
                "variation": f_variation,
                "density": f_density,
                "hardness": f_hardness }

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
  (PF_SLIDER, "f_length", "Fibers Length", 50, (1, 100, 1)),
  (PF_SLIDER, "f_width", "Fibers Width", 50, (1, 100, 1)),
  (PF_SLIDER, "f_orientation", "Fibers Orientation", 0, (-90, 90, 1)),
  (PF_SLIDER, "f_variation", "Fibers Variation", 0, (0, 10, 1)),
  (PF_SLIDER, "f_density", "Fibers Density", 5, (1, 10, 1)),
  (PF_SLIDER, "f_hardness", "Fibers Hardness", 8, (0, 10, 1)),
  (PF_SLIDER, "overlap_gain", "Overlap Gain", 0, (0, 10, 1)),
 ],
 [],
 peignage)

main()
