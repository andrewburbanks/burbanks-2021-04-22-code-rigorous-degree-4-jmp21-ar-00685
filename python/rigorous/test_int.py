import random

import pytest


def random_int():
    a = random.randint(-1000, 1000)
    return a


num_random_tests = 100000


def test_int_o0_notation():
    for k in range(num_random_tests):
        a = random_int()
        b = random_int()
        if k == 0:
            b = int(a)
        if a <= b:
            assert b >= a
            assert not a > b
        elif a >= b:
            assert b <= a
            assert not a < b
        else:
            assert a == b
            assert not a > b
            assert not a < b
            assert not b > a
            assert not b < a


# complete ordered lattice

def test_int_o1_reflexive():
    for k in range(num_random_tests):
        a = random_int()
        assert a <= a


def test_int_o2_transitive():
    for k in range(num_random_tests):
        a = random_int()
        b = random_int()
        c = random_int()
        if a <= b and b <= c:
            assert a <= c


def test_int_o3_antisymmetric():
    for k in range(num_random_tests):
        a = random_int()
        b = random_int()
        if k == 0:
            b = int(a)
        if a <= b:
            if b <= a:
                assert a == b
        elif b <= a:
            if a <= b:
                assert a == b
        else:
            assert a == b


def test_int_o4_total_order():
    for k in range(num_random_tests):
        a = random_int()
        b = random_int()
        assert a <= b or b <= a


def test_int_o5_lattice():
    for k in range(num_random_tests):
        a = random_int()
        b = random_int()
        inf_ab = min(a, b)
        sup_ab = max(a, b)
        assert inf_ab <= a and inf_ab <= b
        assert inf_ab == a or inf_ab == b
        assert a <= sup_ab and b <= sup_ab
        assert a == sup_ab or b == sup_ab


def fail_test_int_o6_conditional_complete_order():
    # how to check this?
    assert False


def fail_test_int_o7_complete_order():
    # how to check this? would have to use -inf, inf as evidence.
    assert False


def fail_test_int_p1_add_associative():
    for k in range(num_random_tests):
        a = random_int()
        b = random_int()
        c = random_int()
        assert (a + b) + c == a + (b + c)


def test_int_p2_add_identity():
    for k in range(num_random_tests):
        a = random_int()
        assert a + 0.0 == a == 0.0 + a


def test_int_p3_add_inverse():
    for k in range(num_random_tests):
        a = random_int()
        assert a + (-a) == 0.0 == (-a) + a


def test_int_p4_add_commutes():
    for k in range(num_random_tests):
        a = random_int()
        b = random_int()
        assert a + b == b + a


def fail_test_int_p5_mul_associative():
    for k in range(num_random_tests):
        a = random_int()
        b = random_int()
        c = random_int()
        assert (a * b) * c == a * (b * c)


def test_int_p6_mul_identity():
    for k in range(num_random_tests):
        a = random_int()
        assert a * 1.0 == a
        assert a == 1.0 * a


def fail_test_int_p7_mul_inverse():
    for k in range(num_random_tests):
        a = random_int()
        if a != 0.0:
            b = 1.0 / a
            assert b == a ** -1
            assert a * b == 1.0, a * b
            assert 1.0 == b * a


def test_int_p8_mul_commutes():
    for k in range(num_random_tests):
        a = random_int()
        b = random_int()
        assert a * b == b * a


def fail_test_int_p9_mul_distributes_add():
    for k in range(num_random_tests):
        a = random_int()
        b = random_int()
        c = random_int()
        assert a * (b + c) == (a * b) + (a * c)


def test_int_p10_trichotomy():
    for k in range(num_random_tests):
        a = random_int()
        if k == 0:
            a = 0.0
        assert a > 0.0 or -a > 0.0 or a == 0.0
        assert not (a > 0.0 and -a > 0.0)
        assert not (a > 0.0 and a == 0.0)
        assert not (-a > 0.0 and a == 0.0)


def test_int_p11_add_closed():
    for k in range(num_random_tests):
        a = random_int()
        b = random_int()
        c = a + b
        assert type(c) is int


def test_int_p12_mul_closed():
    for k in range(num_random_tests):
        a = random_int()
        b = random_int()
        c = a * b
        assert type(c) is int


def fail_test_int_p13_bdd_above_has_lub():
    # how to show this for ints?
    # ints are a finite set
    # we have a maximum operation
    assert False


def fail_test_int_division_is_multiply_reciprocal():
    for k in range(num_random_tests):
        a = random_int()
        b = random_int()
        assert a / b == a * (1.0 / b)


def test_int_ordered_field_add_compatibility():
    for k in range(num_random_tests):
        a = random_int()
        b = random_int()
        c = random_int()
        a, b = sorted((a, b))
        assert a <= b
        assert a + c <= b + c


def test_int_ordered_field_mul_compatibility():
    for k in range(num_random_tests):
        a = random_int()
        b = random_int()
        a, b = abs(a), abs(b)
        assert 0.0 <= a
        assert 0.0 <= b
        assert 0.0 <= a * b
