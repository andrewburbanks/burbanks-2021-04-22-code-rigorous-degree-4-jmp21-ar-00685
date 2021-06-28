from .interval import Interval, RR
from .rectangle import Rectangle

import math

def cover_interval(lo, hi, n):
    """
    Return a list of n overlapping subintervals covering the interval [lo, hi].
    """
    assert lo <= hi
    assert n >= 1
    params_vals = [Interval(k)/Interval(n) for k in range(n+1)]
    sample_vals = [p*Interval(hi) + (Interval(1)-p)*Interval(lo) for p in params_vals]
    x_vals = [Interval.hull([s0, s1]) for s0, s1 in zip(sample_vals[:-1], sample_vals[1:])]
    return x_vals

def cover_graph(f, lo, hi, n):
    """
    Return a list of n overlapping rectangles covering the graph of f on interval [lo, hi].
    """
    assert lo <= hi
    assert n >= 1
    x_vals = cover_interval(lo, hi, n)
    f_vals = [f(x) for x in x_vals]
    z_vals = [Rectangle(x, y) for x, y in zip(x_vals, f_vals)]
    return z_vals

def cover_circle(c, r, m, eps=2**-8):
    """
    Return a list of n=4m overlapping rectangles covering a circle centre c, radius r.
    The parameter eps controls the amount of overlap.
    """
    eps = RR(eps)
    assert m >= 1
    assert eps >= 0
    n = 4*m
    samples = generate_samples(c, r, n, eps=eps)
    grains = combine_samples(samples)
    return grains

def generate_samples(c, r, n, eps=2**-8):
    """
    Generate a list of rectangles, called sample rectangles.
    these are rectangles each guaranteed to contain a point
    on the circle of radius r centered on c (real).
    n+1 rectangles are generated.
    degenerate rectangles along the 4 half-axes are included.
    """
    eps = RR(eps)
    I, R = Interval, Rectangle
    assert n > 0
    assert n%4 == 0
    samples = []
    for k in range(n+1):
        # guarantee that endpoints of monotonic segments are contained
        if k == 0 or k == n:
            x0, y0, x1, y1 = r*(1-eps), 0, r*(1+eps), 0
        elif k == n//4:
            x0, y0, x1, y1 = 0, r*(1-eps), 0, r*(1+eps)
        elif k == n//2:
            x0, y0, x1, y1 = -r*(1-eps), 0, -r*(1+eps), 0
        elif k == n//4*3:
            x0, y0, x1, y1 = 0, -r*(1-eps), 0, -r*(1+eps)
        else:
            theta = 2*RR(math.pi)*k/n
            x0, y0 = r*RR(math.cos(theta))*(1-eps), r*RR(math.sin(theta))*(1-eps)
            x1, y1 = r*RR(math.cos(theta))*(1+eps), r*RR(math.sin(theta))*(1+eps)
        assert I(x0, x0)**2 + I(y0, y0)**2 < I(r, r)**2
        assert I(x1, x1)**2 + I(y1, y1)**2 > I(r, r)**2
        z = R(I(c, c)) + R(I(min(x0, x1), max(x0, x1)),
                           I(min(y0, y1), max(y0, y1)))
        samples.append(z)
    return samples

def combine_samples(samples):
    """
    Combine adjacent sample rectangles to produce a covering of a monotonic curve segment.
    """
    I, R = Interval, Rectangle
    grains = []
    zp = None
    for z in samples:
        if zp is not None:
            grain = R(I(min(z.re.lo, zp.re.lo), max(z.re.hi, zp.re.hi)),
                      I(min(z.im.lo, zp.im.lo), max(z.im.hi, zp.im.hi)))
            grains.append(grain)
        zp = z
    return grains

def sample_line_segment(p0, p1, n):
    """
    Sample n (small) rectangles containing points on the line segment joining p0 to p1.
    The endpoints p0 and p1 are guaranteed to be covered.
    """
    samples = [p0]
    for k in range(1, n):
        a = k/n
        b = Interval(a, a)
        z = p1*b + p0*(Interval(1, 1)-b)
        samples.append(z)
    samples.append(p1)
    return samples

def cover_line_segment(p0, p1, n):
    """
    Return a list of n rectangles covering the line segment joining p0 to p1.
    """
    samples = sample_line_segment(p0, p1, n)
    grains = combine_samples(samples)
    return grains

def cover_polygon(points, n):
    """
    Return a list of rectangles covering a polygon using n points per side.
    """
    augmented_points = points+[points[0]]
    grains = []
    for p0, p1 in zip(augmented_points[:-1], augmented_points[1:]):
        grains = grains + cover_line_segment(p0, p1, n)
    return grains

def rectangles_coords(z_vals):
    """
    Return two separate lists of x, y coordinates for a list of rectangles.
    The output is suitable for the plot command.
    """
    n = len(z_vals)
    x = [list((0.0,)*n) for _ in range(5)]
    y = [list((0.0,)*n) for _ in range(5)]
    for k, z in enumerate(z_vals):
        for j, a in enumerate([z.re.lo, z.re.lo, z.re.hi, z.re.hi, z.re.lo]):
            x[j][k] = a
        for j, a in enumerate([z.im.lo, z.im.hi, z.im.hi, z.im.lo, z.im.lo]):
            y[j][k] = a
    return x, y
