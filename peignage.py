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

# Constants for fiber drawing
MAX_L_VARIATION = .25

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def rotate(self, cos_t, sin_t):
        x = self.x
        y = self.y
        self.x = + x*cos_t + y*sin_t
        self.y = - x*sin_t + y*cos_t

    def translate(self, tx, ty):
        self.x += tx
        self.y += ty

    def is_inside(self, r, wr, hr):
        return (r <= self.x <= wr) and (r <= self.y <= hr)

    def copy(self):
        return Point(self.x, self.y)

class Fiber:
    def __init__(self, point, l, messiness):
        self.l_base = l
        self.m = point
        if messiness == 0:
            self.n = self.m.copy()
            self.n.translate(l, 0)
        else:
            # Displace first point
            dl = random.triangular(0, min(messiness, 5)/10.0, 0)
            dt = random.uniform(-math.pi, math.pi)
            self.m.translate(l*dl*math.cos(dt), l*dl*math.sin(dt))
            # Compute second point
            dl = random.triangular(-MAX_L_VARIATION, MAX_L_VARIATION)
            dt = random.triangular(-math.pi/2.0, math.pi/2.0)
            dt *= messiness**1.6 / 10**1.6
            self.n = self.m.copy()
            self.n.translate(l*(1+dl)*math.cos(dt), l*(1+dl)*math.sin(dt))

    def rotate(self, cos_t, sin_t):
        self.m.rotate(cos_t, sin_t)
        self.n.rotate(cos_t, sin_t)

    def translate(self, tx, ty):
        self.m.translate(tx, ty)
        self.n.translate(tx, ty)

    def trim_to_borders(self, w_img, h_img, r):
        r *= 4          #XXX debug, to be deleted(?)
        wr = w_img - r
        hr = h_img - r
        # both points inside frame
        if self.m.is_inside(r, wr, hr) and self.n.is_inside(r, wr, hr):
            return True
        # both points outside frame & in same external sector
        if ((self.m.x < r and self.n.x < r) or
            (self.m.y < r and self.n.y < r) or
            (self.m.x > wr and self.n.x > wr) or
            (self.m.y > hr and self.n.y > hr)):
            return False
        # case x = constant (...approximately)
        if abs(self.n.x - self.m.x) < 1:
            c = self.m.x
            if self.n.y > self.m.y:
                if self.n.y > hr:
                    self.n.y = hr
                if self.m.y < r:
                    self.m.y = r
            else:
                if self.m.y > hr:
                    self.m.y = hr
                if self.n.y < r:
                    self.n.y = r
            return self.is_long_enough()
        # case y = constant (...approximately)
        if abs(self.n.y - self.m.y) < 1:
            c = self.m.y
            if self.n.x > self.m.x:
                if self.n.x > wr:
                    self.n.x = wr
                if self.m.x < r:
                    self.m.x = r
            else:
                if self.m.x > wr:
                    self.m.x = wr
                if self.n.x < r:
                    self.n.x = r
            return self.is_long_enough()
        # case y = ax + b
        a = (self.n.y - self.m.y) / (self.n.x - self.m.x)
        b = self.m.y - a * self.m.x
        y1 = a*r+b; y2 = a*wr+b; x1 = (r-b)/a; x2 = (hr-b)/a;
        intersections = []
        if r <= y1 <= hr:
            intersections.append((r, y1))
        if r <= y2 <= hr:
            intersections.append((wr, y2))
        if r <= x1 <= wr:
            intersections.append((x1, r))
        if r <= x2 <= wr:
            intersections.append((x2, hr))
        # there should be exactly 0 or 2 intersections
        if not intersections:
            return False
        elif len(intersections) == 2:
            m = Point(intersections[0][0], intersections[0][1])
            n = Point(intersections[1][0], intersections[1][1])
            if self.contains(m):
                if self.contains(n):
                    self.m = m
                    self.n = n
                else:
                    if not self.m.is_inside(r, wr, hr):
                        self.m = m
                    else:
                        self.n = m
            else:
                if not self.m.is_inside(r, wr, hr):
                    self.m = n
                else:
                    self.n = n
            return self.is_long_enough()

    def contains(self, point):
        # test if a point on the fiber line is actually on the fiber segment
        kp = (self.n.y-self.m.y)*(point.y-self.m.y) + \
                (self.n.x-self.m.x)*(point.x-self.m.x)
        kf = self.get_length_2()
        return 0 <= kp <= kf

    def is_long_enough(self):
        return self.get_length() > (1 - MAX_L_VARIATION) * self.l_base

    def get_length_2(self):
        return (self.n.y-self.m.y)**2 + (self.n.x-self.m.x)**2

    def get_length(self):
        return math.sqrt(self.get_length_2())

class Grid:
    def __init__(self, img, fparams):
        # Imported constants
        self.w_img = img.width
        self.h_img = img.height
        self.l = fparams["length"]
        self.d = fparams["width"]
        self.params = fparams
        # Internal frame
        self.r = fparams["width"] / 2.0
        self.w_in = self.w_img - 2*self.r
        self.h_in = self.h_img - 2*self.r
        # Precomputed trigonometrics
        self.cos_t = math.cos(math.radians(fparams["orientation"]))
        self.sin_t = math.sin(math.radians(fparams["orientation"]))
        # Dimensions
        self.w = self.w_in*self.cos_t + self.h_in*abs(self.sin_t) + self.l
        self.h = self.h_in*self.cos_t + self.w_in*abs(self.sin_t)
        # Cells
        self.n_lin = int(self.h / self.d) + 3
        self.points_x0 = [i * self.l / self.n_lin for i in range(self.n_lin)]
        self.n_col = int(self.w / (self.l*1.7)) + 1
        # Translation
        self.tx = self.r - self.l*self.cos_t
        self.ty = self.r + self.l*self.sin_t
        if fparams["orientation"] >= 0:
            self.tx -= self.h_in*self.sin_t*self.cos_t
            self.ty += self.h_in*self.sin_t*self.sin_t
        else:
            self.tx += self.w_in*self.sin_t*self.sin_t
            self.ty += self.w_in*self.sin_t*self.cos_t

        self.reset()

    def reset(self):
        self.points = []
        self.fibers = []
        points_y0 = -self.d * random.random()
        random.shuffle(self.points_x0)
        for i in range(self.n_lin):
            y = points_y0 + i*self.d
            for j in range(self.n_col):
                x = self.points_x0[i] + j*(self.l*1.7)
                self.points.append(Point(x, y))

    def rotate_points(self):
        for point in self.points:
            point.rotate(self.cos_t, self.sin_t)

    def rotate_fibers(self):
        for fiber in self.fibers:
            fiber.rotate(self.cos_t, self.sin_t)

    def translate_points(self):
        for point in self.points:
            point.translate(self.tx, self.ty)

    def translate_fibers(self):
        for fiber in self.fibers:
            fiber.translate(self.tx, self.ty)

    def points_to_fibers(self):
        # Give every point a sibling point to make a fiber
        self.fibers = []
        for point in self.points:
            f = Fiber(point, self.l, self.params["messiness"])
            self.fibers.append(f)

    def trim_to_borders(self):
        # Remove fibers outside borders or too short once trimmed
        self.fibers[:] = [fiber for fiber in self.fibers \
            if fiber.trim_to_borders(self.w_img, self.h_img, self.r) ]

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
    g = Grid(img, fparams)
    g.points_to_fibers()
    g.rotate_fibers()
    g.translate_fibers()
    g.trim_to_borders()

    for fiber in g.fibers:
        layer = random.choice(flayers)
        pdb.gimp_paintbrush_default(layer, 4, [fiber.m.x, fiber.m.y,
                                               fiber.n.x, fiber.n.y])

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
             f_messiness=0,
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
                "messiness": f_messiness,
                "density": f_density,
                "hardness": f_hardness }

    pdb.gimp_image_undo_group_start(img)

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

    # Loop until non-transparent area is large enough
    # XXX use density parameter here

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

    pdb.gimp_image_undo_group_end(img)


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
  (PF_SLIDER, "f_messiness", "Fibers Messiness", 0, (0, 10, 1)),
  (PF_SLIDER, "f_density", "Fibers Density", 5, (1, 10, 1)),
  (PF_SLIDER, "f_hardness", "Fibers Hardness", 8, (0, 10, 1)),
  (PF_SLIDER, "overlap_gain", "Overlap Gain", 0, (0, 10, 1)),
 ],
 [],
 peignage)

main()
