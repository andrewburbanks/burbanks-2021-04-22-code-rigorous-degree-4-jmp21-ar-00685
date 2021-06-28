from .interval import Interval
from .trunc import Trunc
from .ball import Ball
from .function import Domain, Function

import pytest
import random

from .test_trunc import random_float_interval_trunc

num_test_cases = 10000

def test_constant_evaluation():
    b = Ball(Trunc(2, [Interval(1, 1), Interval(0, 0), Interval(0, 0)]), 0.0, 0.0)
    f = Function(Domain(0, 1), b)

    assert f(0) == Interval(1, 1)
    assert f(1) == Interval(1, 1)
    assert f(Interval(0, 0)) == Interval(1, 1)
    assert f(Interval(1, 1)) == Interval(1, 1)
    assert f(Interval(0.0, 0.0)) == Interval(1, 1)
    assert f(Interval(1.0, 1.0)) == Interval(1, 1)

def test_constant_arithmetic():
    b = Ball(Trunc(2, [Interval(1, 1), Interval(0, 0), Interval(0, 0)]), 0.0, 0.0)
    f = Function(Domain(0, 1), b)

    assert f+f == Function(Domain(0, 1), Ball(Trunc(2, [Interval(2, 2), Interval(0, 0), Interval(0, 0)]), 0, 0))
    assert f-f == Function(Domain(0, 1), Ball(Trunc(2, [Interval(0, 0), Interval(0, 0), Interval(0, 0)]), 0, 0))

    assert 2*f == Function(Domain(0, 1), Ball(Trunc(2, [Interval(2, 2), Interval(0, 0), Interval(0, 0)]), 0, 0))
    assert f*2 == Function(Domain(0, 1), Ball(Trunc(2, [Interval(2, 2), Interval(0, 0), Interval(0, 0)]), 0, 0))

    assert f*f == Function(Domain(0, 1), Ball(Trunc(2, [Interval(1, 1), Interval(0, 0), Interval(0, 0)]), 0, 0))

def test_constant_composition():
    b = Ball(Trunc(2, [Interval(1, 1), Interval(0, 0), Interval(0, 0)]), 0.0, 0.0)
    f = Function(Domain(0, 1), b)

    assert f(f) == Function(Domain(0, 1), Ball(Trunc(2, [Interval(1, 1), Interval(0, 0), Interval(0, 0)]), 0, 0))

    assert f.diff_compose(0.5*f) == Function(Domain(0, 1), Ball(Trunc(2, [Interval(0, 0), Interval(0, 0), Interval(0, 0)]), 0, 0))

@pytest.mark.parametrize('domain', [Domain(0, 1), Domain(1, 1), Domain(0, 2.5), Domain(1, 2.5)])
def test_add_constant(domain):
    for k in range(num_test_cases):
        t = random_float_interval_trunc(5)
        b = Ball(t, 0.0, 0.0)
        f = Function(domain, b)
        
        assert f+0 == f
        assert f+0.0 == f
        assert 0+f == f
        assert 0.0+f == f
        
        assert f-0 == f
        assert f-0.0 == f
        assert 0-f == -f
        assert 0.0-f == -f
            
@pytest.mark.parametrize('domain', [Domain(0, 1), Domain(1, 1), Domain(0, 2.5), Domain(1, 2.5)])
def test_mul_constant(domain):
    for k in range(num_test_cases):
        t = random_float_interval_trunc(5)
        b = Ball(t, 0.0, 0.0)
        f = Function(domain, b)
        
        assert f*1 == f
        assert 1*f == f
        assert f*1.0 == f
        assert 1.0*f == f

def test_domain_shift():
    num_local = int(num_test_cases**(1/3)+1)
    for k in range(num_local):
        t = random_float_interval_trunc(5)
        b = Ball(t, 0, 0)
        for j in range(num_local):
            # we will test only for radii that are powers of 2 to allow exact comparison
            c, r = random.randint(-5, 5), random.choice([1, 2, 4, 8])
            f = Function(Domain(c, r), b)

            for i in range(num_local):
                x = random.random()
                assert f(x) == b((x-c)/r)

