#!/usr/bin/python
import math
import random

from gimpfu import pdb


# Constants for fiber drawing
MAX_L_VARIATION = .25
FIBER_INLINE_GAP = 1.5
FRAME_PADDING = .5


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

    def has_transparent_end(self, layer):
        # Detect degraded fibers from a transparent end
        nbr_channels, m = pdb.gimp_drawable_get_pixel(layer, self.m.x, self.m.y)
        nbr_channels, n = pdb.gimp_drawable_get_pixel(layer, self.n.x, self.n.y)
        if nbr_channels < 4:
            return False
        return (m[3] == 0) or (n[3] == 0)

    def get_length_2(self):
        return (self.n.y-self.m.y)**2 + (self.n.x-self.m.x)**2

    def get_length(self):
        # Standard euclidian norm
        return math.sqrt(self.get_length_2())


class Grid:
    """
    A Grid is made of multiple Fibers spread across the frame.
    The fibers can be reset with a call to setup_fibers().
    """
    def __init__(self, img, fparams):
        # Imported constants
        self.w_img = img.width
        self.h_img = img.height
        self.l = fparams["length"]
        self.d = fparams["width"]
        self.params = fparams
        # Internal frame for the fiber points
        self.pad = fparams["width"] / 2.0 * (1 + FRAME_PADDING)
        self.w_in = self.w_img - 2*self.pad
        self.h_in = self.h_img - 2*self.pad
        # Precomputed trigonometrics
        self.cos_t = math.cos(math.radians(fparams["orientation"]))
        self.sin_t = math.sin(math.radians(fparams["orientation"]))
        # Note that once drawn, fibers actually go width/2 further
        overframe = fparams["width"] / 2.0 * fparams["hardness"] / 10.0
        self.max_area = (self.w_in + overframe) * (self.h_in + overframe)
        # Also, there may be two corner areas unreachable to trimmed fibers
        min_l = (1 - MAX_L_VARIATION) * self.l
        min_l_area = min_l * fparams["width"] * fparams["hardness"] / 10.0
        corners_area = self.cos_t * self.sin_t * min_l**2
        corners_area *= (10 - fparams["messiness"]) / 10.0
        corners_area = max(0, corners_area - min_l_area)
        self.max_area -= corners_area
        # Points field dimensions
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

    def decimate(self, x):
        # Drop fibers to lower density
        self.fibers[:] = self.fibers[int((1-x)*len(self.fibers)):]


