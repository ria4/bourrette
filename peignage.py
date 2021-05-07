#!/usr/bin/python
import math
import random
from gimpfu import *


# Hey! It's a debug switch!
DEBUG = False

# Constants for exponential input conversion
L_A = -0.000344125449991
L_B = 0.11153346297
L_C = -7.0189446165
W_A = 3.27285167629e-05
W_B = 0.0594682389829
W_C = -6.96725624648

# Constants for fiber drawing
MAX_L_VARIATION = .25
FIBER_INLINE_GAP = 1.5
FRAME_PADDING = 1.5
LEFTOVERS_THRESHOLD = 200


class Point:
    """
    Simple Point class with basic transformations.
    """
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def rotate(self, cos_t, sin_t):
        # The y axis inversion is accounted for
        x = self.x
        y = self.y
        self.x = + x*cos_t + y*sin_t
        self.y = - x*sin_t + y*cos_t

    def translate(self, tx, ty):
        self.x += tx
        self.y += ty

    def is_inside(self, p, wp, hp):
        # Inside the rectangle (p,p)-(wp,p)-(wp,hp)-(p,hp)
        return (p <= self.x <= wp) and (p <= self.y <= hp)

    def copy(self):
        return Point(self.x, self.y)


class Fiber:
    """
    A Fiber is made of two correlated Points.
    The class builds each instance from one single Point,
    which gets copied then translated with constrained randomness.
    """
    def __init__(self, point, l, messiness):
        self.l_base = l
        self.m = point
        if messiness == 0:
            # Put N exactly l pixels away from M, on the right
            self.n = self.m.copy()
            self.n.translate(l, 0)
        else:
            # Move M around initial position
            dl = random.triangular(0, min(messiness, 5)/10.0, 0)
            dt = random.uniform(-math.pi, math.pi)
            self.m.translate(l*dl*math.cos(dt), l*dl*math.sin(dt))
            # Put N approx. l pixels away, approx. on the right
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

    def trim(self, w_img, h_img, p):
        # If necessary & if possible, trim the fiber to fit inside the frame
        # Return True iff the resulting fiber is deemed to be long enough
        wp = w_img - p
        hp = h_img - p
        # Both points inside frame
        if self.m.is_inside(p, wp, hp) and self.n.is_inside(p, wp, hp):
            # no need to check the length inherited from __init__
            return True
        # Both points outside frame & in same external sector
        if ((self.m.x < p and self.n.x < p) or
            (self.m.y < p and self.n.y < p) or
            (self.m.x > wp and self.n.x > wp) or
            (self.m.y > hp and self.n.y > hp)):
            return False
        # Case x = constant (or close enough)
        if abs(self.n.x - self.m.x) < 1:
            c = self.m.x
            if self.n.y > self.m.y:
                if self.n.y > hp:
                    self.n.y = hp
                if self.m.y < p:
                    self.m.y = p
            else:
                if self.m.y > hp:
                    self.m.y = hp
                if self.n.y < p:
                    self.n.y = p
            return self.is_long_enough()
        # Case y = constant (or close enough)
        if abs(self.n.y - self.m.y) < 1:
            c = self.m.y
            if self.n.x > self.m.x:
                if self.n.x > wp:
                    self.n.x = wp
                if self.m.x < p:
                    self.m.x = p
            else:
                if self.m.x > wp:
                    self.m.x = wp
                if self.n.x < p:
                    self.n.x = p
            return self.is_long_enough()
        # Case y = ax + b, with a =/= 0
        a = (self.n.y - self.m.y) / (self.n.x - self.m.x)
        b = self.m.y - a * self.m.x
        y1 = a*p+b; y2 = a*wp+b; x1 = (p-b)/a; x2 = (hp-b)/a;
        intersections = []
        if p <= y1 <= hp:
            intersections.append((p, y1))
        if p <= y2 <= hp:
            intersections.append((wp, y2))
        if p <= x1 <= wp:
            intersections.append((x1, p))
        if p <= x2 <= wp:
            intersections.append((x2, hp))
        # There should be exactly 0 or 2 intersections
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
                    if not self.m.is_inside(p, wp, hp):
                        self.m = m
                    else:
                        self.n = m
            else:
                if not self.m.is_inside(p, wp, hp):
                    self.m = n
                else:
                    self.n = n
            return self.is_long_enough()

    def contains(self, point):
        # Check whether the given Point is on the Fiber segment
        # NB: This point already OUGHT TO* be on the Fiber extended line
        # *key word as defined in RFC 6919
        kp = (self.n.y-self.m.y)*(point.y-self.m.y) + \
                (self.n.x-self.m.x)*(point.x-self.m.x)
        kf = self.get_length_2()
        return 0 <= kp <= kf

    def is_long_enough(self):
        # Identify fibers shortened beyond 'tasteful' use
        #TODO This is waste. Should we build something else with it?
        return self.get_length() > (1 - MAX_L_VARIATION) * self.l_base

    def get_length_2(self):
        return (self.n.y-self.m.y)**2 + (self.n.x-self.m.x)**2

    def get_length(self):
        # Standard euclidian norm
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
        self.pad = fparams["width"] / 2.0 * FRAME_PADDING
        self.w_in = self.w_img - 2*self.pad
        self.h_in = self.h_img - 2*self.pad
        # Precomputed trigonometrics
        self.cos_t = math.cos(math.radians(fparams["orientation"]))
        self.sin_t = math.sin(math.radians(fparams["orientation"]))
        # Dimensions
        self.w = self.w_in*self.cos_t + self.h_in*abs(self.sin_t) + self.l
        self.h = self.h_in*self.cos_t + self.w_in*abs(self.sin_t)
        # Cells
        self.n_lin = int(self.h / self.d) + 3
        self.points_x0 = [i * self.l*FIBER_INLINE_GAP / self.n_lin \
                for i in range(self.n_lin)]
        self.n_col = int(self.w / (self.l*FIBER_INLINE_GAP)) + 1
        # Translation
        self.tx = self.pad - self.l*self.cos_t
        self.ty = self.pad + self.l*self.sin_t
        if fparams["orientation"] >= 0:
            self.tx -= self.h_in*self.sin_t*self.cos_t
            self.ty += self.h_in*self.sin_t*self.sin_t
        else:
            self.tx += self.w_in*self.sin_t*self.sin_t
            self.ty += self.w_in*self.sin_t*self.cos_t
        # Points & fibers
        self.points = []
        self.fibers = []

    def setup_fibers(self):
        # Call this routine to fill the grid with fresh fibers
        self.create_points()
        self.create_fibers()
        self.rotate_fibers()
        self.translate_fibers()
        self.trim_fibers()
        self.shuffle_fibers()

    def create_points(self):
        # Initial points are spaced evenly on every line,
        # and every line is spaced evenly from its neighbors
        self.points = []
        points_y0 = -self.d * random.random()
        random.shuffle(self.points_x0)
        for i in range(self.n_lin):
            y = points_y0 + i*self.d
            for j in range(self.n_col):
                x = self.points_x0[i] + j*(self.l*FIBER_INLINE_GAP)
                self.points.append(Point(x, y))

    def create_fibers(self):
        # Give every point a sibling point
        self.fibers = []
        for point in self.points:
            f = Fiber(point, self.l, self.params["messiness"])
            self.fibers.append(f)

    def rotate_fibers(self):
        for fiber in self.fibers:
            fiber.rotate(self.cos_t, self.sin_t)

    def translate_fibers(self):
        for fiber in self.fibers:
            fiber.translate(self.tx, self.ty)

    def trim_fibers(self):
        # Remove fibers outside borders or too short once trimmed
        self.fibers[:] = [fiber for fiber in self.fibers \
            if fiber.trim(self.w_img, self.h_img, self.pad) ]

    def shuffle_fibers(self):
        # Mitigate patterns inherited from create_points
        random.shuffle(self.fibers)

    #XXX def decimer


def get_coverage(drawable, area_covered_max):
    # Return the ratio of (alpha-weighted) pixels over the padded frame
    mean, std_dev, median, pixels, count, percentile = \
            pdb.gimp_drawable_histogram(drawable, HISTOGRAM_VALUE, 0, 1)
    return pixels * 1.0 / area_covered_max


def apply_mask_and_merge(img, layer_base, layer_mask, layer_collage, soft_overlap):

    # Create channel based on brightness levels
    # White noise means full effect, black noise means no effect
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
    g = Grid(img, fparams)
    area_covered_max = g.w_in * g.h_in

    # Loop until non-transparent area is large enough
    # XXX use density parameter here

    # Fill grid with proper fibers
    g.setup_fibers()

    # Compute number of fibers per buffer
    n_buffer = len(g.fibers) / n_layers**2
    buffers_end = n_buffer * n_layers**2
    n_leftovers = len(g.fibers) - buffers_end

    # First loop level to prevent layer prioritization
    for i in range(n_layers):

        # Cycle base layers
        base_layers = base_layers[1:] + base_layers[:1]

        # Second loop level to prevent layer prioritization
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
                pdb.gimp_paintbrush_default(layer_mask, 4, coords)

            # Apply layer mask and merge to collage layer
            apply_mask_and_merge(img, layer_base, layer_mask,
                                 layer_collage, fparams["overlap"])

    # Add leftover fibers when it makes a difference
    if (fparams["messiness"] == 0) or (len(g.fibers) < LEFTOVERS_THRESHOLD):

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
            pdb.gimp_paintbrush_default(layer_mask, 4, coords)

            # Apply layer mask and merge to collage layer
            apply_mask_and_merge(img, layer_base, layer_mask,
                                 layer_collage, fparams["overlap"])

    # Clean up
    img.remove_layer(layer_mask)


def peignage(img,
             drawable,
             f_length=50,
             f_width=50,
             f_orientation=0,
             f_messiness=0,
             f_density=6,
             f_hardness=8,
             f_overlap=0):

    # Scale user input exponentially
    f_width_exp = img.height * math.exp(W_A*f_width**2 + W_B*f_width + W_C)
    f_length_exp = img.width * math.exp(L_A*f_length**2 + L_B*f_length + L_C)

    # Store all fiber parameters
    fparams = { "length": f_length_exp,
                "width": f_width_exp,
                "orientation": f_orientation,
                "messiness": f_messiness,
                "density": f_density,
                "hardness": f_hardness,
                "overlap": f_overlap }

    if not DEBUG:
        pdb.gimp_image_undo_group_start(img)

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
        pdb.gimp_image_undo_group_end(img)


register(
 "python-fu-peignage",
 "Mingle layers like silk waste fibers",
 "Mingle layers like silk waste fibers",
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
  (PF_SLIDER, "f_density", "Fibers Density", 6, (0, 10, 1)),
  (PF_SLIDER, "f_hardness", "Fibers Hardness", 8, (0, 10, 1)),
  (PF_SLIDER, "f_overlap", "Soft Overlap", 0, (0, 1, 1)),
  # i don't care much for the bulky PF_BOOL button
 ],
 [],
 peignage)

main()
