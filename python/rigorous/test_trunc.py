import pytest
import random

from .trunc import Trunc
from .test_interval import random_float_interval
from .interval import Interval, RR

max_degree = 10
num_cases = 100


def random_float():
    a = random.randint(-1000, 1000)
    b = random.random()
    return RR(a) + RR(b)


def random_int():
    a = random.randint(-1000, 1000)
    return a


def random_float_trunc(n):
    coeffs = [random_float() for k in range(n + 1)]
    return Trunc(n, coeffs)


def random_float_interval_trunc(n):
    coeffs = [random_float_interval() for k in range(n + 1)]
    return Trunc(n, coeffs)


def random_int_trunc(n):
    coeffs = [random_int() for k in range(n + 1)]
    return Trunc(n, coeffs)


@pytest.mark.parametrize('n', range(11))
def test_make_zero(n):
    p = Trunc(n, [0] * (n + 1))
    assert repr(p) == 'Trunc(%d, [%s])' % (n, ', '.join(['0'] * (n + 1)))


@pytest.mark.parametrize('n', range(11))
def test_make_fewer_coeffs(n):
    assert Trunc(1, [1, 0, 0]) == Trunc(1, [1, 0]) == Trunc(1, [1])
    assert Trunc(2, [1, 0, 0]) == Trunc(2, [1, 0]) == Trunc(2, [1])
    assert Trunc(2, [0, 1, 0]) == Trunc(2, [0, 1])
    for m in range(1, n + 1):
        for k in range(num_cases):
            coeffs = [random_float() for j in range(m + 1)]
            p = Trunc(n, coeffs)
            q = Trunc(n, coeffs + [0.0] * (n - m))
            assert p == q


@pytest.mark.parametrize('n', range(11))
def test_additive_identity(n):
    zero = Trunc(n, [0] * (n + 1))
    for k in range(num_cases):
        p = random_float_trunc(n)
        assert p + zero == p
        assert zero + p == p


@pytest.mark.parametrize('n', range(1, 11))
def test_compose_one_left_identity(n):
    for k in range(num_cases):
        p = random_float_trunc(n)
        assert 1 * p == p


@pytest.mark.parametrize('n', range(1, 11))
def test_compose_one_right_identity(n):
    for k in range(num_cases):
        p = random_float_trunc(n)
        assert p * 1 == p


@pytest.mark.parametrize('n', range(1, 11))
def test_compose_left_identity(n):
    x = Trunc(n, [0, 1] + [0] * (n - 1))
    for k in range(num_cases):
        p = random_float_trunc(n)
        assert x(p) == p, x(p)


@pytest.mark.parametrize('n', range(1, 11))
def test_compose_right_identity(n):
    x = Trunc(n, [0, 1] + [0] * (n - 1))
    for k in range(num_cases):
        p = random_float_trunc(n)
        assert p(x) == p


@pytest.mark.parametrize('n', range(1, 11))
def test_eval_zero_gives_constant_term(n):
    for k in range(num_cases):
        p = random_float_trunc(n)
        assert p(0) == p[0]


@pytest.mark.parametrize('n', range(1, 11))
def test_eval_one_sums_coeffs(n):
    for k in range(num_cases):
        p = random_float_trunc(n)
        assert p(1) == sum(p.coeffs)




@pytest.mark.parametrize('n', range(2, 11))
def test_interval_coeffs(n):
    p = Trunc(n, [Interval(0, 0), Interval(1, 1)] + [Interval(0, 0) for k in range(n - 1)])
    for k in range(n + 1):
        assert p ** k == Trunc.basis_element(n, k, Interval(1, 1))


@pytest.mark.parametrize('n', range(1, 11))
def test_basis_elements(n):
    one = Interval(1, 1)
    zero = one * 0
    for k in range(n + 1):
        e_k = Trunc.basis_element(n, k, one)
        for j in range(n + 1):
            if j == k:
                assert e_k[j] == one
            else:
                assert e_k[j] == zero


@pytest.mark.parametrize('n', range(1, 11))
def test_basis_elements_mul(n):
    one = 1
    zero = one * 0
    for j in range(n + 1):
        e_j = Trunc.basis_element(n, j, one)
        for k in range(n + 1):
            e_k = Trunc.basis_element(n, k, one)
            if j + k <= n:
                assert e_j * e_k == Trunc.basis_element(n, j + k, one)
            else:
                assert e_j * e_k == Trunc.constant(n, zero)


@pytest.mark.parametrize('n', range(1, 11))
def test_trunc_is_coeff_dot_basis_elements(n):
    for k in range(num_cases):
        p = random_float_trunc(n)
        zero = p * 0
        total = zero
        for j in range(n + 1):
            total = total + p[j] * Trunc.basis_element(n, j, p[0] ** 0)
        assert total == p


# ITERATION

@pytest.mark.parametrize('n', range(1, 11))
def test_trunc_iterable_list(n):
    for k in range(num_cases):
        p = random_float_trunc(n)
        coeffs = list(p)


@pytest.mark.parametrize('n', range(1, 11))
def test_trunc_iterable_trunc_equals_trunc_list_coeffs(n):
    for k in range(num_cases):
        p = random_float_trunc(n)
        assert p == Trunc(n, list(p))


@pytest.mark.parametrize('n', range(1, 11))
def test_trunc_iterable_for(n):
    for k in range(num_cases):
        p = random_float_trunc(n)
        for c in p:
            pass


@pytest.mark.parametrize('n', range(1, 11))
def test_trunc_iterable_for_enumerate(n):
    for k in range(num_cases):
        p = random_float_trunc(n)
        for k, c in enumerate(p):
            assert p[k] == c


# DERIVATIVE

def test_trunc_derivative_examples_constant():
    assert Trunc(0, [1]).diff() == Trunc(0, [0])
    assert Trunc(0, [1.0]).diff() == Trunc(0, [0.0])


def test_trunc_derivative_examples_constant_interval():
    assert Trunc(0, [Interval(1, 1)]).diff() == Trunc(0, [Interval(0, 0)])
    assert Trunc(0, [Interval(1.0, 1.0)]).diff() == Trunc(0, [Interval(0.0, 0.0)])


def test_trunc_derivative_examples():
    assert Trunc(5, [0, 1, 2, 3, 4, 5]).diff() == Trunc(5, [1, 4, 9, 16, 25])


@pytest.mark.parametrize('deg', range(1, 11))
def test_trunc_derivative_pure_powers(deg):
    assert Trunc(deg, [0] * deg + [1]).diff() == Trunc(deg, [0] * (deg - 1) + [deg])


@pytest.mark.parametrize('n', range(1, 11))
def test_trunc_derivative_linear_exact(n):
    # pure powers and linearity should be enough to verify polynomial derivative.
    for k in range(num_cases):
        p = random_int_trunc(n)
        q = random_int_trunc(n)
        a = random_int()
        b = random_int()
        assert (p + q).diff() == p.diff() + q.diff()
        assert (a * p + b * q).diff() == a * p.diff() + b * q.diff()
