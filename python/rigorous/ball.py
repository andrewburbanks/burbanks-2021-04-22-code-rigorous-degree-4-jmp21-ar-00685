from .interval import RR, Interval, Pub
from .rectangle import Rectangle
from .poly import Poly
from .trunc import Trunc

def bound_quasi_power(k0, g_norm):
    assert g_norm < Interval(1.0, 1.0)
    k = k0 + 1
    latest = k * g_norm ** (k - 1)
    best = latest
    while latest >= best:
        best = latest
        k = k + 1
        latest = k * g_norm ** (k - 1)
    return best


class Ball:

    def __init__(self, P, H, G):
        assert isinstance(P, Trunc)
        assert all(isinstance(a, Interval) for a in P)
        self.P = P
        self.H = Pub(H) # Ball.__init__ ensures interval [0, hi]
        self.G = Pub(G) # Ball.__init__ ensures interval [0, hi]

    def __repr__(self):
        return 'Ball(P=%r, H=%r, G=%r)' % (self.P, self.H, self.G)

    @property
    def truncation_degree(self):
        return self.P.truncation_degree

    def __eq__(self, other):
        # for comparing balls, not the functions contained within them
        assert isinstance(other, Ball)
        return self.P == other.P and self.H == other.H and self.G == other.G

    def __pos__(self):
        P = +self.P
        H = self.H
        G = self.G
        return Ball(P, H, G)

    def __neg__(self):
        P = -self.P
        H = self.H
        G = self.G
        return Ball(P, H, G)

    def add_bi(self, other):
        assert isinstance(other, Interval)
        other = Ball.constant(self.truncation_degree, other)
        return self.add_bb(other)

    def add_bb(self, other):
        assert isinstance(other, Ball)
        P = self.P + other.P
        H = self.H + other.H
        G = self.G + other.G
        return Ball(P, H, G)

    def __add__(self, other):
        if isinstance(other, Ball):
            return self.add_bb(other)
        if isinstance(other, Interval):
            return self.add_bi(other)
        if isinstance(other, (RR, float, int)):
            return self.add_bi(Interval(other, other))
        assert False

    def sub_bi(self, other):
        assert isinstance(other, Interval)
        other = Ball.constant(self.truncation_degree, other)
        return self.sub_bb(other)

    def sub_bb(self, other):
        assert isinstance(other, Ball)
        P = self.P - other.P
        H = self.H + other.H
        G = self.G + other.G
        return Ball(P, H, G)

    def __sub__(self, other):
        if isinstance(other, Ball):
            return self.sub_bb(other)
        if isinstance(other, Interval):
            return self.sub_bi(other)
        if isinstance(other, (RR, float, int)):
            return self.sub_bi(Interval(other, other))
        if isinstance(other, Trunc):
            return self.sub_bb(Ball(other, 0, 0))
        assert False, f'{type(other)}'

    def __mul__(self, other):
        if isinstance(other, int) and other == 0:
            return Ball(self.P * other, self.H * other, self.G * other)
        if isinstance(other, (Interval, RR, float, int)):
            return Ball(self.P * other, self.H * abs(other), self.G * abs(other))
        if isinstance(other, Trunc):
            return self*Ball(other, 0.0, 0.0)
        P = self.P * other.P
        H = (Trunc.bound_high_order_product(self.P, other.P)
             + abs(self.P) * other.H + self.H * abs(other.P)
             + self.H * other.H
             + self.H * other.G + self.G * other.H)
        G = (abs(self.P) * other.G + self.G * abs(other.P)
             + self.G * other.G)
        return Ball(P, H, G)

    def __rmul__(self, other):
        return self * other

    def __pow__(self, p):
        assert isinstance(p, int)
        assert p >= 0
        zero = Ball(self.P * 0, self.H * 0, self.G * 0)
        one = Ball(self.P ** 0, self.H * 0, self.G * 0)
        prod = one
        power = self
        while p > 0:
            if (p % 2) == 1:
                prod = prod * power
            power = power * power
            p = p // 2
        return prod

    def __call__(self, other):
        # this can be made much more efficient
        # no need to create all three balls

        if isinstance(other, Ball):
            return self.comp_bb(other)
        if isinstance(other, Trunc):
            assert self.truncation_degree == other.truncation_degree
            return self(Ball(other, 0.0, 0.0))
        if isinstance(other, Poly):
            if other.degree <= self.truncation_degree:
                fP, H, G = other, 0.0, 0.0
            else:
                fP, fH = other.split(self.truncation_degree)
                H = abs(fH)
            return self(Ball(Trunc(self.truncation_degree, list(fP)), H, G))
        if isinstance(other, (Rectangle, Interval)):
            return self.comp_bj(other)
        if isinstance(other, (RR, float, int)):
            return self(Interval(other, other))
        raise NotImplementedError

    def comp_bj(self, other):
        assert isinstance(other, (Rectangle, Interval)), repr(other)
        if self.H == self.G == Interval(0, 0):
            pass
        else:
            assert abs(other) in Interval(0, 1)
        N = self.P.truncation_degree
        g_norm = abs(other) ** (N + 1)
        r = self.H * g_norm + self.G
        return self.P(other) + Interval(-(r.hi), r.hi)

    def comp_bj_slow(self, other):
        assert isinstance(other, Interval)
        b = Ball.constant(self.truncation_degree, other)
        ans = self.comp_bb(b)
        assert ans.H == Interval(0.0, 0.0)
        for k in range(1, self.truncation_degree + 1):
            assert ans.P[k] == Interval(0.0, 0.0)
        return Interval((ans.P[0] - ans.G).lo, (ans.P[0] + ans.G).hi)

    def comp_bb(self, other):
        assert isinstance(other, Ball)
        #print('Composition F o G with ||G||_1 in ', abs(other))
        # the following test is not needed if H,G are both Pub(0.0).
        if self.H == self.G == Interval(0, 0):
            pass
        else:
            assert abs(other) in Interval(0, 1), 'Composition F o G undefined as ||G||_1 >= 1'
        P_zero = self.P * 0
        ball_1 = self.P(other)
        N = self.P.truncation_degree
        g_norm = abs(other) ** (N + 1)
        g_star_norm = Ball.abs_star(other) ** (N + 1)
        ball_2 = Ball(P_zero, self.H * g_star_norm, self.H * (g_norm - g_star_norm))
        ball_3 = Ball(P_zero, self.H * 0, self.G)
        return ball_1 + ball_2 + ball_3

    def __abs__(self):
        x = abs(self.P) + Interval(0, self.H.hi) + Interval(-self.G.hi, self.G.hi)
        return Interval(max(0 * x.lo, x.lo), x.hi)

    def abs_r(self):
        hi = abs(self).hi
        return Interval(hi, hi)

    def abs_star(self):
        x = Trunc.abs_star(self.P) + Interval(0, self.H.hi) + Interval(-self.G.hi, self.G.hi)
        return Interval(max(0 * x.lo, x.lo), x.hi)

    def diff_compose(self, other):
        #print('Ball.diff_compose')
        # this can be made much more efficient
        # no need to create all three balls
        if isinstance(other, Trunc):
            return self.diff_compose(Ball(other, 0.0, 0.0))

        assert isinstance(other, Ball), f'Ball.diff_compose attempt with arg {type(other)}'
        if self.G == self.H == Interval(0, 0):
            pass
        else:
            assert abs(other) in Interval(0, 1)
        P_zero = self.P * 0

        ball_1 = self.P.diff()(other)
        #assert ball_1.H == Interval(0, 0)  # not true in general
        #assert ball_1.G == Interval(0, 0)  # not true in general

        N = self.P.truncation_degree
        g_norm = other.abs_r()
        #print(other)
        #print(other.abs_r())
        if self.G == self.H == Interval(0, 0):
            norm_sup_H = norm_sup_G = 0
        else:
            norm_sup_H = bound_quasi_power(self.truncation_degree, g_norm)
            norm_sup_G = bound_quasi_power(0, g_norm)
        #print('Bounds on quasi powers', norm_sup_H, norm_sup_G)
        # some sacrifice is made here! absorbing into general term!
        # an alternative is to absord into both high-order term and
        # the final polynomial coefficient.
        ball_2 = Ball(P_zero, self.H * 0, self.H * norm_sup_H)
        ball_3 = Ball(P_zero, self.H * 0, self.G * norm_sup_G)
        return ball_1 + ball_2 + ball_3

    def __contains__(self, other):
        assert isinstance(other, Poly)
        if isinstance(other[0], Interval):
            raise NotImplementedError
            # Containment only implemented for numerical coeffs
        P, H = other.split(self.P.truncation_degree)
        if all(a in I for a, I in zip(P, self.P)):
            return abs(H) in self.H + self.G
        else:
            return abs(H) in self.H + (self.G - self.bound_norm_distance_to_poly(P))

    def bound_norm_distance_to_poly(self, P):
        assert isinstance(P, Poly)
        assert not isinstance(P[0], Interval)
        zero = self.P[0] * 0
        total = zero
        for a, I in zip(P, self.P):
            if a < I.lo:
                norm_increment_to_shift_a_into_I = Interval(I.lo, I.lo) - Interval(a, a)
            elif a > I.hi:
                norm_increment_to_shift_a_into_I = Interval(a, a) - Interval(I.hi, I.hi)
            else:
                assert a in I
                norm_increment_to_shift_a_into_I = zero
            assert norm_increment_to_shift_a_into_I >= Interval(0, 0)
            total = total + norm_increment_to_shift_a_into_I
        return total

    @staticmethod
    def basis_element(n, k, one):
        assert 0 <= k <= n + 1
        if k <= n:
            return Ball(Trunc.basis_element(n, k, one), 0.0, 0.0)
        else:
            return Ball(Trunc.constant(n, one * 0), 1.0, 0.0)

    @staticmethod
    def constant(n, a):
        return Ball.basis_element(n, 0, a ** 0) * a

