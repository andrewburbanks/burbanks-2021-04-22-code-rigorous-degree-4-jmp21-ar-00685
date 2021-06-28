import random

import pytest

def random_float():
    a = random.randint(-1000, 1000)
    b = random.random()
    return float(a) + b


num_random_tests = 100000


def test_float_o0_notation():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        if k == 0:
            b = float(a)
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


"""
The reals, R, when extended with -oo and +oo to Rhat,
form a complete lattice under <=.

The finite floats, F, when extended with -inf and +inf to Fhat,
form a screen of the above lattice.  (For this, we associate
-inf with -oo and +inf with +oo.)

Recall that a screen S of a set M is a subset such that
(1) every element of M has lower and upper bounds in S,
(2a) the set of upper bounds for each elt a of M has a least elt in S,
(2b) the set of lower bounds for each elt a of M has a greatest elt in S,
(3) the screen is symmetric iff M has unary minus op such that a in S implies -a in S.

The above are true for floats because <= for floats is identical to <= for reals
restricted to floats, and because the sign-bit allows exact unary minus.

The upshot of this is that a screen forms a complete sublattice with the same least
and greatest element (at least, when we identify infinite floats with the corresponding
extended reals).  The idea is that this enables some of the structure of, e.g., the reals,
to carry-over to the floats.
 
A map from the complete lattice Rhat onto the screen Fhat is called a Rounding if
(R1) forall a in S: round(a) == a, (ADB: i.e. round is a projection onto S). 
(R2) a<=b implies round(a)<=round(b), (ADB: i.e., round is order preserving).

A rounding is directed down, resp. directed up, if
(R3) forall a in M: round(a) <= a, resp., round(a) >= a.

A rounding is antisymmetric if exists unary minus on M and
(R4) round(-a) == -round(a),
NOTE: standard roundings to nearest or towards or away from zero are antisymmetric.

Monotone directed downward (resp. upward) roundings of a complete lattice onto a screen are unique.

The point of all of this is that we want to carry properties and ops over to floats.

Consider a complete lattice {M,<=} and S a screen of M, op_M a binop on M,
then a rounding round can be used to approximate op_M by op_S defined by
(RG) a op_S b := round(a op_M b).
If exists unary minus on M, and S is symmetric screen (wrt unary minus)
then (R1,2,4) and (RG) give a Semimorphism.

Such semimorphisms where the roundings happen to be antisymmetric are very well-suited
to transferring properties from M to S.

e.g., Ordered Field and Ordered Vector Space axioms are preserved.

If a, b are adjacent floats in Fhat with a <= x <= b and x in Rhat,
and round is a rounding, then (R1) gives a <= round(x) <= b,
so there are never any other floats between x and round(x).

Similary, for the operations {+,-,*,/}, we have
a <= x +_Rhat y <= b then by (R1,2) and (RG) a <= x +_Fhat y <= b,
provided that x +_Fhat y := round(x +_Rhat y).  IEEE754 should ensure this!

The properties of floats can even be defined to be those of Rhat that
are preserved by semimorphism, via antisymmetric rounding.

NOTE:
Extending the above to intervals needs care with the meaning of <=
which really means subseteq for closed (non-extended but possibly unbounded)
real intervals and
under which these intervals form a complete lattice again.

With this meaning, closed (non-extended) float intervals form a screen of the complete lattice.

With this meaning, monotone directed upwards rounding means superseteq
rounding, i.e., lower bounds round down, upper bounds round up.  This is
the key to understanding interval containment.

Thus, axioms for reals/floats have <= meaning less-than-or-equal,
those for real intervals / float intervals have <= meaning subset-or-equal.

IEEE754 ensures that directed rounding modes on floats (up and down) can be
combined correctly to achieve the result of monotone upward (wrt subseteq)
rounding on the intervals, and hence the most important property of inclusion.

Could check
(IR1)
(IR2)
(IR3)
(IR4)

(OD1)
(OD2)
(OD3)
(OD4)
(OD5)
(IRG)

There are some technicalities around exception-free calculus by avoiding using
extended reals for closed intervals.
"""


# complete ordered lattice

def test_float_o1_reflexive():
    for k in range(num_random_tests):
        a = random_float()
        assert a <= a


def test_float_o2_transitive():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        c = random_float()
        if a <= b and b <= c:
            assert a <= c


def test_float_o3_antisymmetric():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        if k == 0:
            b = float(a)
        if a <= b:
            if b <= a:
                assert a == b
        elif b <= a:
            if a <= b:
                assert a == b
        else:
            assert a == b


def test_float_o4_total_order():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        assert a <= b or b <= a


def test_float_o5_lattice():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        inf_ab = min(a, b)
        sup_ab = max(a, b)
        assert inf_ab <= a and inf_ab <= b
        assert inf_ab == a or inf_ab == b
        assert a <= sup_ab and b <= sup_ab
        assert a == sup_ab or b == sup_ab


def fail_test_float_o6_conditional_complete_order():
    # how to check this?
    assert False


def fail_test_float_o7_complete_order():
    # how to check this? would have to use -inf, inf as evidence.
    assert False


def fail_test_float_p1_add_associative():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        c = random_float()
        assert (a + b) + c == a + (b + c)


def test_float_p2_add_identity():
    for k in range(num_random_tests):
        a = random_float()
        assert a + 0.0 == a == 0.0 + a


def test_float_p3_add_inverse():
    for k in range(num_random_tests):
        a = random_float()
        assert a + (-a) == 0.0 == (-a) + a


def test_float_p4_add_commutes():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        assert a + b == b + a


def fail_test_float_p5_mul_associative():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        c = random_float()
        assert (a * b) * c == a * (b * c)


def test_float_p6_mul_identity():
    for k in range(num_random_tests):
        a = random_float()
        assert a * 1.0 == a
        assert a == 1.0 * a


def fail_test_float_p7_mul_inverse():
    for k in range(num_random_tests):
        a = random_float()
        if a != 0.0:
            b = 1.0 / a
            assert b == a ** -1
            assert a * b == 1.0, a * b
            assert 1.0 == b * a


def test_float_p8_mul_commutes():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        assert a * b == b * a


def fail_test_float_p9_mul_distributes_add():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        c = random_float()
        assert a * (b + c) == (a * b) + (a * c)


def test_float_p10_trichotomy():
    for k in range(num_random_tests):
        a = random_float()
        if k == 0:
            a = 0.0
        assert a > 0.0 or -a > 0.0 or a == 0.0
        assert not (a > 0.0 and -a > 0.0)
        assert not (a > 0.0 and a == 0.0)
        assert not (-a > 0.0 and a == 0.0)


def test_float_p11_add_closed():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        c = a + b
        assert type(c) is float


def test_float_p12_mul_closed():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        c = a * b
        assert type(c) is float


def fail_test_float_p13_bdd_above_has_lub():
    # how to show this for floats?
    # floats are a finite set
    # we have a maximum operation
    assert False


def fail_test_float_division_is_multiply_reciprocal():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        assert a / b == a * (1.0 / b)


def test_float_ordered_field_add_compatibility():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        c = random_float()
        a, b = sorted((a, b))
        assert a <= b
        assert a + c <= b + c


def test_float_ordered_field_mul_compatibility():
    for k in range(num_random_tests):
        a = random_float()
        b = random_float()
        a, b = abs(a), abs(b)
        assert 0.0 <= a
        assert 0.0 <= b
        assert 0.0 <= a * b
