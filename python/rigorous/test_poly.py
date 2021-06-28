import pytest
import random

from .poly import Poly, trim
from .test_interval import random_float_interval
from .interval import Interval, RR

max_degree = 10
num_cases = 100

"""
Polynomials in X with coefficients from K is denoted K[X].
Assumed that K forms a field?

R is a ring iff: (2 ops on same set with distrib)
Abelian group under +, with 0 as additive identity:
(1) + associative
(2) + identity 0
(3) + inverses
(4) + commutative
Monoid under *:
(5) * associative
(6) * identity
Distributative:
(7) * distributes over +.

F is a field iff:
Abelian group under +, with 0 as additive identity:
(1) + associative
(2) + commutative
(3) + identity 0
(4) + inverses
Nonzero elements abelian group under *:
(5) * associative
(6) * commutative
(7) * identity (distinct from identity for +)
(8) * inverses except for 0 [ensures trivial (single elt) ring is not a field]
Distributative:
(9) * distributes over +.

The reals form a field. The floats do not. The ints do not. The rationals do.
We could test using rationals as coeffs?

Given a field for K, polynomials K[X] form a commutative algebra.
Using a commutative ring instead of a field delivers a polynomial ring as above.

A module is an abelian group (say under +) together with a ring (scalars) such
that scalar multiplication is defined. For scalars r, s: (two sets)
(1) s*(x+y) = s*x+s*y,
(2) (r+s)*x = r*x+s*x,
(3) (rs)*x = r*(s*x)
(4) 1_R*x = x where 1_R is the ring (scalar) identity.
[In fact, the above is a left R-module.]

An algebra is a module M together with a multiplication compatible with the scalar
multiplication. 

Graded algebra so that grades of a and b determine the grade of ab and we can decompose.
Relevant for multivariate polynomials.
"""


# CREATE RANDOM TEST CASES

def random_float():
    a = random.randint(-10, 10)
    b = random.random()
    return RR(a) + RR(b)


def random_int():
    a = random.randint(-10, 10)
    return a


def random_float_poly(n):
    coeffs = [random_float() for k in range(n + 1)]
    return Poly(coeffs)


def random_int_poly(n):
    coeffs = [random_int() for k in range(n + 1)]
    return Poly(coeffs)


def random_float_interval_poly(n):
    coeffs = [random_float_interval() for k in range(n + 1)]
    return Poly(coeffs)



# BASIC HELPER FUNCTIONS

def test_poly_trim_examples():
    assert trim([0, 1, 2, 3, 4]) == [0, 1, 2, 3, 4]
    assert trim([0, 1, 2, 3, 0]) == [0, 1, 2, 3]
    assert trim([0, 1, 2, 0, 0]) == [0, 1, 2]
    assert trim([0, 1, 0, 0, 0]) == [0, 1]
    assert trim([0, 0, 0, 0, 0]) == [0]
    assert trim([0, 1, 0, 3, 4]) == [0, 1, 0, 3, 4]
    assert trim([0, 1, 0, 3, 0]) == [0, 1, 0, 3]
    assert trim([0, 1, 0, 0, 0]) == [0, 1]
    assert trim([0, 1, 0, 0, 0]) == [0, 1]
    assert trim([0, 0, 0, 0, 0]) == [0]


def test_poly_trim_interval_examples():
    # This test ensures that trim performs comparisons with
    # zero that also work for intervals.
    i0 = Interval(0, 0)
    i1 = Interval(1, 1)
    i2 = Interval(2, 2)
    i3 = Interval(3, 3)
    i4 = Interval(4, 4)
    assert trim([i0, i1, i2, i3, i4]) == [i0, i1, i2, i3, i4]
    assert trim([i0, i1, i2, i3, i0]) == [i0, i1, i2, i3]
    assert trim([i0, i1, i2, i0, i0]) == [i0, i1, i2]
    assert trim([i0, i1, i0, i0, i0]) == [i0, i1]
    assert trim([i0, i0, i0, i0, i0]) == [i0]
    assert trim([i0, i1, i0, i3, i4]) == [i0, i1, i0, i3, i4]
    assert trim([i0, i1, i0, i3, i0]) == [i0, i1, i0, i3]
    assert trim([i0, i1, i0, i0, i0]) == [i0, i1]
    assert trim([i0, i1, i0, i0, i0]) == [i0, i1]
    assert trim([i0, i0, i0, i0, i0]) == [i0]


# ZERO

@pytest.mark.parametrize('n', range(11))
def test_poly_zero_make(n):
    p = Poly([0] * (n + 1))
    assert repr(p) == 'Poly([0])'
    # assert repr(p) == 'Poly1d(%d, [%s])' % (n, ', '.join(['0'] * (n + 1)))


def test_poly_zero_degree():
    """
    This is a design choice.  Strictly-speaking, -inf would be a good degree.
    """
    z = Poly([0])
    assert z.degree == 0


def test_poly_zero_mul_zero():
    z = Poly([0])
    assert z * z == z


def test_poly_zero_add_zero():
    z = Poly([0])
    assert z + z == z


def test_poly_zero_sub_zero():
    z = Poly([0])
    assert z - z == z


def test_poly_zero_pow_int_zero():
    z = Poly([0])
    o = Poly([1])
    assert z ** 0 == o


# def test_FAIL_poly_zero_compares_zero():
#    zero = Poly1d([0])
#    assert zero == Poly1d([0])
#    assert zero == Poly1d([0.0])
#    assert zero == Poly1d([Interval(0, 0)])
#    assert zero == Poly1d([Interval(0.0, 0.0)])


# RING AXIOMS

@pytest.mark.parametrize('n', range(11))
def test_poly_ring_1_add_associative_exact(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        q = random_int_poly(n)
        r = random_int_poly(n)
        assert (p + q) + r == p + (q + r)


@pytest.mark.parametrize('n', range(11))
def test_poly_ring_2_add_identity(n):
    zero = Poly([0])
    for k in range(num_cases):
        p = random_float_poly(n)
        assert p + zero == p
        assert zero + p == p


@pytest.mark.parametrize('n', range(11))
def test_poly_ring_3_add_inverse_exact(n):
    zero = Poly([0])
    for k in range(num_cases):
        p = random_int_poly(n)
        q = -p
        assert p + q == zero
        assert q + p == zero


@pytest.mark.parametrize('n', range(11))
def test_poly_ring_4_add_commutes(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        q = random_float_poly(n)
        assert p + q == q + p


@pytest.mark.parametrize('n', range(11))
def test_poly_ring_5_mul_associative_exact(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        q = random_int_poly(n)
        r = random_int_poly(n)
        assert (p * q) * r == p * (q * r)


@pytest.mark.parametrize('n', range(11))
def test_poly_ring_6_mul_identity(n):
    one = Poly([1])
    for k in range(num_cases):
        p = random_float_poly(n)
        assert p * one == p
        assert one * p == p


@pytest.mark.parametrize('n', range(11))
def test_poly_ring_7_distributative_exact(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        q = random_int_poly(n)
        r = random_int_poly(n)
        assert p * (q + r) == (p * q) + (p * r)
        assert (p + q) * r == (p * r) + (q * r)


# UNARY OPERATORS

@pytest.mark.parametrize('n', range(11))
def test_poly_pos(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        assert +p == p


@pytest.mark.parametrize('n', range(11))
def test_poly_pos_new(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        assert +p is not p


@pytest.mark.parametrize('n', range(11))
def test_poly_neg_involution(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        assert --p == p


@pytest.mark.parametrize('n', range(11))
def test_poly_neg_new(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        assert --p is not p


# SUBTRACTION

@pytest.mark.parametrize('n', range(11))
def test_poly_sub_self(n):
    zero = Poly([0])
    for k in range(num_cases):
        p = random_int_poly(n)
        assert p - p == zero


@pytest.mark.parametrize('n', range(11))
def test_poly_sub_additive_inverse_exact(n):
    zero = Poly([0])
    for k in range(num_cases):
        p = random_int_poly(n)
        q = random_int_poly(n)
        assert p - q == p + (-q)


@pytest.mark.parametrize('n', range(11))
def test_poly_subtraction_additive_inverse_inexact_coeffs(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        q = random_float_poly(n)
        assert p - q == p + (-q)


# @pytest.mark.parametrize('n', range(11))
# def test_mul_FAIL_commutes_inexact_Coeffs(n):
#    """
#    Float polynomial multiplication does not commute here,
#    due to failures of associativity for float addition,
#    and failure of distributivity.
#    """
#    for k in range(num_cases):
#        p = random_float_poly(n)
#        q = random_float_poly(n)
#        assert p * q == q * p


@pytest.mark.parametrize('n', range(11))
def test_poly_ring_7_mul_commutes_exact(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        q = random_int_poly(n)
        assert p * q == q * p


# MODULE AXIOMS
# add abelian group:
# test_poly_module_1_add_associative already tested for ring_1
# test_poly_module_2_add_identity already tested for ring_2
# test_poly_module_3_add_inverses already tested for ring_3
# test_poly_module_4_add_commutative already tested for ring_4
# scalar mul:
# 5-8 match vector space axioms
#
# A few thoughts:
# If K is a field then a K-vector space is the same as a K-module.

@pytest.mark.parametrize('n', range(1, 11))
def test_poly_module_5_scalar_mul_distributes_vector_add_exact(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        q = random_int_poly(n)
        a = random_int()
        assert a * (p + q) == (a * p) + (a * q)


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_module_6_scalar_mul_distributes_scalar_add(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        a = random_int()
        b = random_int()
        assert (a + b) * p == (a * p) + (b * p)


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_module_7_scalar_mul_compatible_with_field_mul_exact(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        a = random_int()
        b = random_int()
        assert a * (b * p) == (a * b)*p


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_module_8_scalar_mul_left_identity(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        assert RR(1) * p == p


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_module_8a_scalar_mul_right_identity(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        assert p * RR(1) == p


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_module_scalar_mul_commutes(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        a = random_float()
        assert p * a == a * p


# VECTOR SPACE AXIOMS with scalars from field F
# test_poly_vector_1_add_associative already tested for ring_1
# test_poly_vector_2_add_identity already tested for ring_2
# test_poly_vector_3_add_inverses already tested for ring_3
# test_poly_vector_4_add_commutative already tested for ring_4
# test_poly_vector_6_scalar_mul_identity already tested for module
# the remaining axioms are to do with scalar multiplication only
# these would fail for float coeffs, although would succeed within tolerance

@pytest.mark.parametrize('n', range(1, 11))
def test_poly_vector_5_scalar_mul_compatible_with_field_mul_exact(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        a = random_int()
        b = random_int()
        assert a * (b * p) == (a * b) * p


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_vector_7_scalar_mul_distributes_vector_add_exact(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        q = random_int_poly(n)
        a = random_int()
        assert a * (p + q) == (a * p) + (a * q)


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_vector_8_scalar_mul_distributes_scalar_add(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        a = random_int()
        b = random_int()
        assert (a + b) * p == (a * p) + (b * p)


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_vector_scalar_mul_commutes(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        a = random_float()
        assert a * p == p * a


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_vector_scalar_mul_zero_gives_zero_vector(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        assert RR(0.0) * p == Poly([0.0])


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_vector_scalar_mul_negone_gives_neg(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        assert RR(-1.0) * p == -p


# COMPOSITION (monoid)
# - associative
# - identity
# - no inverses in general
# - does not commute in general

@pytest.mark.parametrize('n', range(1, 6))
def test_poly_compose_associative(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        q = random_int_poly(n)
        r = random_int_poly(n)
        assert p(q)(r) == p(q(r))


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_compose_left_identity(n):
    x = Poly([0, 1])
    for k in range(num_cases):
        p = random_float_poly(n)
        assert x(p) == p, x(p)


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_compose_right_identity(n):
    x = Poly([0, 1])
    for k in range(num_cases):
        p = random_float_poly(n)
        assert p(x) == p


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_compose_left_constant(n):
    for k in range(num_cases):
        a = random_float()
        p = Poly([a])
        q = random_float_poly(n)
        assert p(q) == p


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_compose_right_constant_exact(n):
    for k in range(num_cases):
        a = random_int()
        p = Poly([a])
        q = random_int_poly(n)
        composed = q(p)
        evaluated = q(a)
        assert composed == Poly([evaluated])


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_compose_left_linear_exact(n):
    for k in range(num_cases):
        p = random_int_poly(n)
        q = random_int_poly(n)
        r = random_int_poly(n)
        assert (p + q)(r) == p(r) + q(r)


# EVALUATION

@pytest.mark.parametrize('n', range(1, 11))
def test_poly_eval_zero_gives_constant_term(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        assert p(0) == p[0]


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_eval_one_sums_coeffs(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        assert p(1) == sum(p.coeffs)


def test_poly_eval_constant_gives_constant_term():
    for k in range(num_cases):
        a = random_float()
        p = Poly([a])
        x = random_float()
        assert p(x) == a


# VECTOR SPACE BASIS

@pytest.mark.parametrize('n', range(1, 11))
def test_poly_basis_elements_mul(n):
    one = 1
    zero = one * 0
    for j in range(n + 1):
        e_j = Poly.basis_element(n, j, one)
        for k in range(n + 1):
            e_k = Poly.basis_element(n, k, one)
            assert e_j * e_k == Poly.basis_element(n + n, j + k, one)


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_is_coeff_dot_basis_elements(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        zero = p * RR(0)
        total = zero
        for j in range(n + 1):
            total = total + p[j] * Poly.basis_element(n, j, p[0] ** 0)
        assert total == p


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_basis_elements_pow(n):
    one = 1
    zero = one * 0
    x = Poly([zero, one])
    for j in range(n + 1):
        e_j = Poly.basis_element(n, j, one)
        assert x ** j == e_j


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_basis_elements_make(n):
    one = Interval(1, 1)
    zero = one * 0
    for k in range(n + 1):
        e_k = Poly.basis_element(n, k, one)
        assert e_k.degree == k
        for j in range(k + 1):
            if j == k:
                assert e_k[j] == one
            else:
                assert e_k[j] == zero


@pytest.mark.parametrize('n', range(2, 11))
def test_poly_basis_interval(n):
    p = Poly([Interval(0, 0), Interval(1, 1)] + [Interval(0, 0) for k in range(n - 1)])
    for k in range(n + 1):
        lhs = p ** k
        rhs = Poly.basis_element(n, k, Interval(1, 1))
        assert lhs == rhs


# TRUNCATION

@pytest.mark.parametrize('n', range(1, 11))
def test_poly_truncation_examples(n):
    for k in range(num_cases):
        coeffs = [random_float() for _ in range(n + 1)]
        p = Poly(coeffs)
        for j in range(k * 2):
            d = min(j, p.degree)
            assert p.truncate(j) == Poly(coeffs[:d + 1])


def test_poly_truncation_explicit():
    assert Poly([0, 1, 2]).truncate(5) == Poly([0, 1, 2])
    assert Poly([0, 1, 2]).truncate(4) == Poly([0, 1, 2])
    assert Poly([0, 1, 2]).truncate(3) == Poly([0, 1, 2])
    assert Poly([0, 1, 2]).truncate(2) == Poly([0, 1, 2])
    assert Poly([0, 1, 2]).truncate(1) == Poly([0, 1])
    assert Poly([0, 1, 2]).truncate(0) == Poly([0])


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_truncated_product(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        q = random_float_poly(n)
        pq = p * q
        for j in range(2 * k + 1):
            assert Poly.truncated_product(p, q, j) == pq.truncate(j)


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_split_identity(n):
    for k in range(num_cases):
        coeffs = [random_float() for _ in range(n + 1)]
        p = Poly(coeffs)
        for j in range(k * 2):
            d = min(j, p.degree)
            p_lo, p_hi = p.split(j)
            assert p_lo.degree <= j
            assert p == p_lo + p_hi


# DEGREE

@pytest.mark.parametrize('n', range(1, 11))
@pytest.mark.parametrize('m', range(1, 11))
def test_poly_mul_degree(n, m):
    for k in range(num_cases):
        p = random_float_poly(n)
        q = random_float_poly(m)
        pq = p * q
        assert pq.degree <= p.degree + q.degree


@pytest.mark.parametrize('n', range(1, 11))
@pytest.mark.parametrize('m', range(1, 11))
def test_poly_mul_degree(n, m):
    for k in range(num_cases):
        p = random_float_poly(n)
        q = random_float_poly(m)
        pq = p * q
        if p * RR(0) != p and q * RR(0) != q:
            assert pq.degree == p.degree + q.degree


@pytest.mark.parametrize('n', range(1, 11))
@pytest.mark.parametrize('m', range(1, 11))
def test_poly_add_degree(n, m):
    for k in range(num_cases):
        p = random_float_poly(n)
        q = random_float_poly(m)
        assert (p + q).degree <= max(p.degree, q.degree)


# CONSTANT POLYNOMIALS

import operator

unary_operators = [operator.pos,
                   operator.neg]

binary_operators = [operator.add,
                    operator.sub,
                    operator.mul]


@pytest.mark.parametrize('op', unary_operators)
def test_poly_constant_unary_ops(op):
    for k in range(num_cases):
        a = random_float()
        p = Poly([a])
        assert op(p) == Poly([op(a)])


@pytest.mark.parametrize('op', binary_operators)
def test_poly_constant_binary_ops(op):
    for k in range(num_cases):
        a = random_float()
        b = random_float()
        p = Poly([a])
        q = Poly([b])
        assert op(p, q) == Poly([op(a, b)])


@pytest.mark.parametrize('n', range(11))
def test_poly_constant_mul_compatible_scalar_mul(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        a = random_float()
        assert a * p == Poly([a]) * p
        assert p * a == p * Poly([a])


# ITERATION

@pytest.mark.parametrize('n', range(1, 11))
def test_poly_iterable_list(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        coeffs = list(p)


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_iterable_poly_equals_poly_list_coeffs(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        assert p == Poly(list(p))


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_iterable_for(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        for c in p:
            pass


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_iterable_for_enumerate(n):
    for k in range(num_cases):
        p = random_float_poly(n)
        for k, c in enumerate(p):
            assert p[k] == c


def test_poly_derivative_examples_constant():
    assert Poly([1]).diff() == Poly([0])
    assert Poly([1.0]).diff() == Poly([0.0])


def test_poly_derivative_examples_constant_interval():
    assert Poly([Interval(1, 1)]).diff() == Poly([Interval(0, 0)])
    assert Poly([Interval(1.0, 1.0)]).diff() == Poly([Interval(0.0, 0.0)])


def test_poly_derivative_examples():
    assert Poly([0, 1, 2, 3, 4, 5]).diff() == Poly([1, 4, 9, 16, 25])


@pytest.mark.parametrize('deg', range(1, 11))
def test_poly_derivative_pure_powers(deg):
    assert Poly([0] * deg + [1]).diff() == Poly([0] * (deg - 1) + [deg])


@pytest.mark.parametrize('n', range(1, 11))
def test_poly_derivative_linear_exact(n):
    # pure powers and linearity should be enough to verify polynomial derivative.
    for k in range(num_cases):
        p = random_int_poly(n)
        q = random_int_poly(n)
        a = random_int()
        b = random_int()
        assert (p + q).diff() == p.diff() + q.diff()
        assert (a * p + b * q).diff() == a * p.diff() + b * q.diff()

