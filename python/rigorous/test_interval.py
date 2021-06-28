# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 09:55:48 2018

@author: Dr Andrew Burbanks
"""

import pytest

from .interval import Interval, RR

import random


def random_float():
    a = random.randint(-1000, 1000)
    b = random.random()
    return RR(a) + RR(b)


def random_float_interval():
    a = random_float()
    b = random_float()
    a, b = sorted((a, b))
    return Interval(a, b)


num_random_tests = 10000
num_samples_per_interval = 5

"""
NOTE: integer arithmetic is exact in Python, so should we use ints to test
basic algebraic properties first, before using floats?

NOTE: float addition is commutative, but not associative!
"""


def test_interval_random():
    for k in range(num_random_tests):
        x = random_float_interval()
        assert isinstance(x, type(Interval(0, 0)))


# COMPARATORS

def test_lt_true_examples():
    # pyinterval does not implement these correctly
    assert Interval(-1, 1) < Interval(2, 3)


def test_lt_false_examples():
    # pyinterval does not implement these correctly
    assert not Interval(-1, 1) < Interval(1, 2)
    assert not Interval(-1, 1) < Interval(0, 1)
    assert not Interval(-1, 1) < Interval(-1, 0)
    assert not Interval(-1, 1) < Interval(-2, -1)
    assert not Interval(-1, 1) < Interval(-3, -2)


def test_gt_false_examples():
    assert not Interval(-1, 1) > Interval(2, 3)
    assert not Interval(-1, 1) > Interval(1, 2)
    assert not Interval(-1, 1) > Interval(0, 1)
    assert not Interval(-1, 1) > Interval(-1, 0)
    assert not Interval(-1, 1) > Interval(-2, -1)


def test_gt_true_examples():
    assert Interval(-1, 1) > Interval(-3, -2)


def test_eq_false_examples():
    # pyinterval does not implement these correctly
    assert not Interval(-1, 1) == Interval(2, 3)
    assert not Interval(-1, 1) == Interval(1, 2)
    assert not Interval(-1, 1) == Interval(0, 1)
    assert not Interval(-1, 1) == Interval(-1, 0)
    assert not Interval(-1, 1) == Interval(-2, -1)


def test_eq_true_examples():
    for k in range(num_random_tests):
        x = random_float_interval()
        assert x == x


# CREATION

def test_interval_create_zero():
    x = Interval(0, 0)
    assert x.lo == 0
    assert x.hi == 0


def test_interval_bounds_stored():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        a, b = sorted((a, b))
        x = Interval(a, b)
        assert x.lo == a
        assert x.hi == b


def test_interval_equals_self():
    for k in range(num_random_tests):
        x = random_float_interval()
        assert x == x


def test_interval_notequals_self_plus_nonzero():
    one = Interval(1, 1)
    for k in range(num_random_tests):
        x = random_float_interval()
        assert x != x + one


def test_interval_addition_commutes():
    for k in range(num_random_tests):
        x = random_float_interval()
        y = random_float_interval()
        assert x + y == y + x


def test_interval_addition_associative_should_fail():
    with pytest.raises(AssertionError):
        for k in range(num_random_tests):
            x = random_float_interval()
            y = random_float_interval()
            z = random_float_interval()
            assert (x + y) + z == x + (y + z)


def test_interval_addition_identity():
    zero = Interval(0.0, 0.0)
    for k in range(num_random_tests):
        x = random_float_interval()
        assert x + zero == zero + x
        assert x + zero == x  # fails for pessimistic bounds


# def test_interval_addition_checks_type():
#    for k in range(num_random_tests):
#        p = random_interval()
#        with pytest.raises(Exception):
#            p+0

def test_interval_pos_exists():
    for k in range(num_random_tests):
        x = random_float_interval()
        +x


def test_interval_neg_exists():
    for k in range(num_random_tests):
        x = random_float_interval()
        -x


def test_interval_addition_negative_commutes():
    for k in range(num_random_tests):
        x = random_float_interval()
        assert x + (-x) == (-x) + x


def test_interval_additive_inverse_fails():
    """This should fail due to the dependency problem!
    It should result in a symmetric interval, though."""
    zero = Interval(0, 0)
    for k in range(num_random_tests):
        x = random_float_interval()
        if x.hi != x.lo:
            assert x + (-x) != zero
            assert (-x) + x != zero


def test_interval_additive_inverse_special():
    zero = Interval(0, 0)
    for k in range(num_random_tests):
        a = random_float()
        x = Interval(a, a)
        b = x + (-x)
        c = (-x) + x
        assert b == c
        assert b == zero
        assert c == zero
        assert b.lo == -b.hi
        assert c.lo == -c.hi

def test_interval_subtraction_neg_compatible():
    for k in range(num_random_tests):
        x = random_float_interval()
        y = random_float_interval()
        assert x - y == -(y - x)


def test_interval_subtraction_of_neg():
    for k in range(num_random_tests):
        x = random_float_interval()
        y = random_float_interval()
        assert x - (-y) == x + y


def test_interval_additive_inverse_contains_symmetric_span():
    for k in range(num_random_tests):
        x = random_float_interval()
        span = x.hi - x.lo
        assert Interval(-span, span) in x + (-x)


def test_interval_subtraction_adds_negative():
    zero = Interval(0, 0)
    for k in range(num_random_tests):
        x = random_float_interval()
        y = random_float_interval()
        assert x - y == x + (-y)


def test_interval_contains_bounds():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        a, b = sorted((a, b))
        x = Interval(a, b)
        assert a in x
        assert b in x


def test_interval_multiply_zero():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        a, b = sorted((a, b))
        x = Interval(a, b)
        assert a * 0 == b * 0 == 0
        assert x * 0 == Interval(0, 0), (a, b)


def test_interval_power_zero():
    """pyinterval has an issue with raising to
    power zero.  If the interval being raised
    contains zero, then the resulting interval is
    [0, 1] instead of [1, 1]."""
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        a, b = sorted((a, b))
        x = Interval(a, b)
        assert a ** 0 == 1
        assert b ** 0 == 1
        assert x ** 0 == Interval(1, 1), (a, b)


def test_interval_power_one():
    for k in range(num_random_tests):
        x = random_float_interval()
        assert x ** 1 == x


def pow_mul(x, n):
    assert n > 0
    p = x
    for k in range(2, n + 1):
        p = p * x
    return p


@pytest.mark.parametrize('p', range(2, 11))
def test_interval_power_examples(p):
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        a, b = sorted((a, b))
        x = Interval(a, b)
        c, d = sorted((a ** p, b ** p))
        assert Interval(c, d) in x ** p


# now tests that are repeated for each operator

import itertools
import operator


def sample_interval(x, n):
    lo, hi = x.lo, x.hi
    for k in range(n):
        p = RR(k) / RR(n - 1)
        a = p * hi + (1 - p) * lo
        if k == 0:
            assert a == lo
        if k == n - 1:
            assert a == hi
        yield a


unary_operators = [operator.pos,
                   operator.neg]

binary_operators = [operator.add,
                    operator.sub,
                    operator.mul,
                    operator.truediv]


@pytest.mark.parametrize('op', unary_operators)
def test_containment_unary(op):
    for k in range(num_random_tests):
        x = random_float_interval()
        a_vals = sample_interval(x, num_samples_per_interval)
        for a in a_vals:
            assert op(a) in op(x)


@pytest.mark.parametrize('op', binary_operators)
def test_containment_binary(op):
    for k in range(num_random_tests):
        x = random_float_interval()
        y = random_float_interval()
        if op == operator.truediv and 0 in y:
            continue
        a_vals = sample_interval(x, num_samples_per_interval)
        b_vals = sample_interval(y, num_samples_per_interval)
        for (a, b) in itertools.product(a_vals, b_vals):
            assert op(a, b) in op(x, y)


@pytest.mark.parametrize('op', unary_operators)
def test_interval_degenerate_acts_like_scalar_unary(op):
    for k in range(num_random_tests):
        x = random_float()
        a = Interval(x, x)
        assert Interval(op(x), op(x)) in op(a)


@pytest.mark.parametrize('op', binary_operators)
def test_interval_degenerate_acts_like_scalar_binary(op):
    for k in range(num_random_tests):
        x = random_float()
        y = random_float()
        a = Interval(x, x)
        b = Interval(y, y)
        # if op == operator.truediv and y == 0:
        #    continue
        assert Interval(op(x, y), op(x, y)) in op(a, b)
