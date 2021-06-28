import pytest
import random

from .ball import Ball
from .poly import Poly
from .trunc import Trunc
from .interval import Interval, RR

from .test_interval import random_float_interval
from .test_poly import random_float_poly, random_int_poly, random_float_interval_poly
from .test_trunc import random_float_interval_trunc

# we intend to test using polys and zero bounds

num_test_cases = 10


def split_poly(p, n):
    assert len(p) >= 1
    p_lo = p.truncate(n)
    p_hi = p - p_lo
    return p_lo, p_hi


def random_float_interval_ball(n):
    p = random_float_interval_trunc(n)
    h = RR(random.random())
    g = RR(random.random())
    b = Ball(p, h, g)
    return b


def diff(seq):
    elts = iter(seq)
    a = next(elts)
    for b in elts:
        yield b - a
        a = b


def random_float_vector_bounded_norm(n, bound=1):
    # n elements needed
    # chosen as differences between n+1 sorted randoms on [0, 1)
    coeffs = list(diff(sorted(RR(random.random()) * bound for k in range(n + 1))))
    return coeffs


def random_float_interval_ball_bounded_norm(n):
    P1, H, G = random_float_vector_bounded_norm(3)
    P = Trunc(n, [Interval(a, a) for a in random_float_vector_bounded_norm(n + 1, P1)])
    b = Ball(P, H, G)
    return b


def test_ball_role_of_constant_term_in_composition():
    # we choose examples in which coeffs and bounds are representable and exactly summable
    # we choose an extreme example in which g has constant term 1.
    f = Ball(Trunc(2, [Interval(1, 1), Interval(2, 2), Interval(3, 3)]), 0.5, 0.25)
    g = Ball(Trunc(2, [Interval(1, 1), Interval(0, 0), Interval(0, 0)]), 0.0, 0.0)
    h = f(g)
    assert h == Ball(Trunc(2, [Interval(6, 6), Interval(0, 0), Interval(0, 0)]), 0.0, 0.75)


def test_ball_poly_vs_ball_composition_example():
    f = Poly([1, 2, 3])
    g = Poly([0.5, 0.125, 0.25])
    h = f(g)
    assert h == Poly([2.75, 0.625, 1.296875, 0.1875, 0.1875])
    #
    bf_tight = Ball(Trunc(1, [Interval(1, 1), Interval(2, 2)]), 3, 0)
    bg_tight = Ball(Trunc(1, [Interval(0.5, 0.5), Interval(0.125, 0.125)]), 0.25, 0)
    bh = bf_tight(bg_tight)
    assert f in bf_tight
    assert g in bg_tight
    assert h in bh
    assert bh == Ball(Trunc(1, [Interval(2, 2), Interval(0.25, 0.25)]), 0.921875, 2.25)
    #
    bh_tight = Ball(Trunc(1, [Interval(2.75, 2.75), Interval(0.625, 0.625)]), 1.296875 + 0.1875 + 0.1875, 0)
    assert h in bh_tight


def test_ball_role_of_high_order_terms_in_composition():
    f = Ball(Trunc(2, [Interval(1, 1), Interval(2, 2), Interval(3, 3)]), 0.5, 0.25)
    g = Ball(Trunc(2, [Interval(0, 0), Interval(0, 0), Interval(0, 0)]), 1.0, 0.0)
    h = f(g)
    assert h == Ball(Trunc(2, [Interval(1, 1), Interval(0, 0), Interval(0, 0)]), 5.5, 0.75)


def test_ball_role_of_general_terms_in_composition():
    f = Ball(Trunc(2, [Interval(1, 1), Interval(2, 2), Interval(3, 3)]), 0.5, 0.25)
    g = Ball(Trunc(2, [Interval(0, 0), Interval(0, 0), Interval(0, 0)]), 0.0, 1.0)
    h = f(g)
    assert h == Ball(Trunc(2, [Interval(1, 1), Interval(0, 0), Interval(0, 0)]), 0.5, 5.75)


def test_ball_identity_in_composition():
    f = Ball(Trunc(2, [Interval(1, 1), Interval(2, 2), Interval(3, 3)]), 0.5, 0.25)
    g = Ball(Trunc(2, [Interval(0, 0), Interval(1, 1), Interval(0, 0)]), 0.0, 0.0)
    h = f(g)
    assert h == Ball(Trunc(2, [Interval(1, 1), Interval(2, 2), Interval(3, 3)]), 0.5, 0.25)


def test_ball_quadratic_in_composition():
    f = Ball(Trunc(2, [Interval(1, 1), Interval(2, 2), Interval(3, 3)]), 0.5, 0.25)
    g = Ball(Trunc(2, [Interval(0, 0), Interval(0, 0), Interval(1, 1)]), 0.0, 0.0)
    h = f(g)
    assert h == Ball(Trunc(2, [Interval(1, 1), Interval(0, 0), Interval(2, 2)]), 3.5, 0.25)


def test_ball_contains_poly_lo_P_examples():
    # need to test some false examples too
    i0 = Interval(0, 0)
    i1 = Interval(1, 1)
    assert Poly([1, 0, 0, 0]) in Ball(Trunc(1, [i1, i0]), 0, 0)
    assert Poly([0, 1, 0, 0]) in Ball(Trunc(1, [i0, i1]), 0, 0)


def test_ball_contains_poly_hi_H_examples():
    i0 = Interval(0, 0)
    i1 = Interval(1, 1)
    assert Poly([0, 0, 1, 0]) in Ball(Trunc(1, [i0, i0]), 1, 0)
    assert Poly([0, 0, 0, 1]) in Ball(Trunc(1, [i0, i0]), 1, 0)


def test_ball_contains_poly_lo_G_examples():
    i0 = Interval(0, 0)
    i1 = Interval(1, 1)
    assert Poly([1, 0, 0, 0]) in Ball(Trunc(1, [i0, i0]), 0, 1)
    assert Poly([0, 1, 0, 0]) in Ball(Trunc(1, [i0, i0]), 0, 1)


def test_ball_contains_poly_hi_G_examples():
    i0 = Interval(0, 0)
    i1 = Interval(1, 1)
    assert Poly([0, 0, 1, 0]) in Ball(Trunc(1, [i0, i0]), 0, 1)
    assert Poly([0, 0, 0, 1]) in Ball(Trunc(1, [i0, i0]), 0, 1)


# DERIVATIVE COMPOSITION

def test_ball_poly_vs_ball_diff_compose_example():
    f = Poly([1, 2, 3, 4])
    g = Poly([0.5, 0.125, 0.25, 0.0625])
    df = f.diff()
    h = df(g)
    assert h == Poly([8.0, 2.25, 4.6875, 1.875, 0.9375, 0.375, 0.046875])
    #
    bf_tight = Ball(Trunc(1, [Interval(1, 1), Interval(2, 2)]), 3 + 4, 0)
    bg_tight = Ball(Trunc(1, [Interval(0.5, 0.5), Interval(0.125, 0.125)]), 0.25 + 0.0625, 0)
    bh = bf_tight.diff_compose(bg_tight)
    print(bh)
    assert f in bf_tight
    assert g in bg_tight
    assert h in bh
    # these bounds will be pretty terrible, with everything in the general term...
    assert bh == Ball(Trunc(1, [Interval(2.0, 2.0), Interval(0.0, 0.0)]), 0.0, RR('42.53898945130751320375939175'))
    bh_tight = Ball(Trunc(1, [Interval(8.0, 8.0), Interval(2.25, 2.25)]), 7.921875, 0)
    assert h in bh_tight


from itertools import product


@pytest.mark.parametrize('n', range(1, 11))
def test_ball_contains_extrema_of_polynomial_boundary(n):
    for k in range(num_test_cases):
        p = random_float_interval_poly(n)
        t = Trunc(len(p) - 1, list(p))
        b = Ball(t, 0, 0)
        lo_bounds = [c.lo for c in p]
        hi_bounds = [c.hi for c in p]
        p_lo = Poly(lo_bounds)
        p_hi = Poly(hi_bounds)
        assert p_lo in b
        assert p_hi in b


@pytest.mark.parametrize('n', range(1, 11))
def test_ball_extrema_of_polynomial_have_distance_zero_to_boundary(n):
    for k in range(num_test_cases):
        p = random_float_interval_poly(n)
        t = Trunc(len(p) - 1, list(p))
        b = Ball(t, 0, 0)
        lo_bounds = [c.lo for c in p]
        hi_bounds = [c.hi for c in p]
        p_lo = Poly(lo_bounds)
        p_hi = Poly(hi_bounds)
        assert b.bound_norm_distance_to_poly(p_lo) == Interval(0, 0)
        assert b.bound_norm_distance_to_poly(p_hi) == Interval(0, 0)


def test_ball_bound_distance_to_poly_example():
    p = Poly([1, 2, 3, 4, 5, 6])
    t = Trunc(3, [Interval(-1, 1), Interval(0, 0), Interval(5, 6), Interval(-3, 5)])
    b = Ball(t, 0, 0)
    assert b.bound_norm_distance_to_poly(p) == Interval(4, 4)
    assert p not in Ball(t, 0, 0)
    assert p not in Ball(t, 5 + 6, 4 - 1)
    assert p not in Ball(t, 5 + 6 - 1, 4)
    assert p in Ball(t, 5 + 6, 4)
    assert p not in Ball(t, 0, 4 + 5 + 6 - 1)
    assert p in Ball(t, 0, 4 + 5 + 6)


@pytest.mark.parametrize('n', range(1, 11))
def test_ball_call_identity(n):
    x = Ball(Trunc(n, [Interval(0, 0), Interval(1, 1)] + [Interval(0, 0)] * (n - 1)), 0, 0)
    for k in range(num_test_cases):
        b = random_float_interval_ball_bounded_norm(n)
        assert x(b) == b
        assert b(x) == b


@pytest.mark.parametrize('n', range(1, 11))
def test_ball_zero_power_zero(n):
    zero = Ball(Trunc(n, [Interval(0, 0)] * (n + 1)), 0, 0)
    one = Ball(Trunc(n, [Interval(1, 1)] + [Interval(0, 0)] * n), 0, 0)
    assert zero ** 0 == one


@pytest.mark.parametrize('n', range(1, 11))
def test_ball_add_identity(n):
    zero = Ball(Trunc(n, [Interval(0, 0)] * (n + 1)), 0, 0)
    for k in range(num_test_cases):
        b = random_float_interval_ball(n)
        assert b + zero == b
        assert zero + b == b


@pytest.mark.parametrize('n', range(1, 11))
def test_ball_mul_absorb(n):
    zero = Ball(Trunc(n, [Interval(0, 0)] * (n + 1)), 0, 0)
    for k in range(num_test_cases):
        b = random_float_interval_ball(n)
        assert zero * b == zero


@pytest.mark.parametrize('n', range(1, 11))
def test_ball_mul_identity(n):
    one = Ball(Trunc(n, [Interval(1, 1)] + [Interval(0, 0)] * n), 0, 0)
    for k in range(num_test_cases):
        b = random_float_interval_ball(n)
        assert b * one == b, f'{b*one}\n{b}'
        assert one * b == b, f'{one*b}\n{b}'


@pytest.mark.parametrize('n', range(1, 11))
def test_ball_init_example(n):
    i1 = Interval(1, 1)
    i0 = i1 * 0
    p = Trunc(n, [i1] + [i0] * n)
    b = Ball(p, 0, 0)
    assert b.truncation_degree == n
    assert b == b


@pytest.mark.parametrize('n', range(1, 11))
def test_ball_random_example(n):
    for k in range(num_test_cases):
        b = random_float_interval_ball(n)


@pytest.mark.parametrize('n', range(0, 11))
def test_ball_call_containment(n):
    for k in range(num_test_cases):
        f = random_float_poly(n)
        # but this should refuse due to norm of g exceeding 1?
        g = Poly(random_float_vector_bounded_norm(n + 1))
        assert abs(g) < 1
        for j in range(1, n + 1):
            f_P, f_H = split_poly(f, j)
            g_P, g_H = split_poly(g, j)
            u_H = abs(f_H)
            v_H = abs(g_H)
            f_T = Trunc(j, [Interval(a, a) for a in f_P])
            g_T = Trunc(j, [Interval(a, a) for a in g_P])
            F = Ball(f_T, u_H, 0)
            G = Ball(g_T, v_H, 0)
            assert f(g) in F(G)


import operator

unary_operators = [operator.pos, operator.neg]
binary_operators = [operator.add, operator.sub, operator.mul]


@pytest.mark.parametrize('n', range(1, 11))
@pytest.mark.parametrize('op', unary_operators)
def test_ball_unary_containment(n, op):
    for k in range(num_test_cases):
        f = random_float_poly(n)
        for j in range(1, n + 1):
            f_P, f_H = split_poly(f, j)
            u_H = abs(f_H)
            f_T = Trunc(j, [Interval(a, a) for a in f_P])
            i0 = Interval(0, 0)
            F = Ball(f_T, u_H, 0)
            assert op(f) in op(F)


@pytest.mark.parametrize('n', range(1, 11))
@pytest.mark.parametrize('op', binary_operators)
def test_ball_binary_containment(n, op):
    for k in range(num_test_cases):
        f = random_int_poly(n)
        g = random_int_poly(n)
        for j in range(1, n + 1):
            f_P, f_H = split_poly(f, j)
            g_P, g_H = split_poly(g, j)
            u_H = abs(f_H)
            v_H = abs(g_H)
            f_T = Trunc(j, [Interval(a, a) for a in f_P])
            g_T = Trunc(j, [Interval(a, a) for a in g_P])
            F = Ball(f_T, u_H, 0)
            G = Ball(g_T, v_H, 0)
            assert op(f, g) in op(F, G)


# we could use general polys and check containment op ops within balls of truncated polys
# this way, we would have access to the high-order coefficients.
# NOTE: we would need to use an exact type (e.g., int or int interval or rational interval) for poly
# in order for all the binary tests to succeed.  This is because intervals resulting from poly
# calculations with increments need not be tight and might not, therefore, leave enough
# norm contribution for successful containment.
