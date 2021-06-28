from itertools import islice
from .interval import Interval, RR


def sumup(seq):
    elts = iter(seq)
    total = next(elts)
    for x in elts:
        total = total + x
    return total


class Trunc:

    # this class could be improved in efficiency
    # we could allow shorter coefficient lists than the truncation degree
    # using zip would allow padding by coeff[0]*0 with zip_longest
    # this would need some amendments to multiplication, but the unit tests
    # would help to ensure that the code didn't break, plus we can keep
    # the slow version for comparison.

    def __init__(self, truncation_degree, coeffs):
        assert truncation_degree >= 0
        assert 1 <= len(coeffs)  # <= truncation_degree + 1
        self.truncation_degree = truncation_degree
        self.coeffs = list(coeffs[:min(len(coeffs), truncation_degree + 1)])
        for k in range(len(self.coeffs), truncation_degree + 1):
            zero = self.coeffs[0] * 0
            self.coeffs.append(zero)
        assert len(self.coeffs) == truncation_degree + 1

    def __repr__(self):
        return 'Trunc(%r, %r)' % (self.truncation_degree, self.coeffs)

    def __len__(self):
        return len(self.coeffs)

    def __getitem__(self, k):
        return self.coeffs[k]

    def __eq__(self, other):
        return (self.truncation_degree == other.truncation_degree
                and all(a == b for a, b in zip(self.coeffs, other.coeffs)))

    def __pos__(self):
        return Trunc(self.truncation_degree, [+a for a in self.coeffs])

    def __neg__(self):
        return Trunc(self.truncation_degree, [-a for a in self.coeffs])

    def __add__(self, other):
        if isinstance(other, (Interval, RR, float, int)):
            return Trunc(self.truncation_degree,
                         [self.coeffs[0] + other] + list(self.coeffs[1:]))
        assert isinstance(other, Trunc), f'{type(other)}'
        assert other.truncation_degree == self.truncation_degree
        return Trunc(self.truncation_degree, [a + b for a, b in zip(self.coeffs, other.coeffs)])

    def __sub__(self, other):
        if isinstance(other, (Interval, RR, float, int)):
            return Trunc(self.truncation_degree,
                         [self.coeffs[0] - other] + list(self.coeffs[1:]))
        assert isinstance(other, Trunc)
        assert other.truncation_degree == self.truncation_degree
        return Trunc(self.truncation_degree, [a - b for a, b in zip(self.coeffs, other.coeffs)])

    def __mul__(self, other):
        if isinstance(other, int) and other == 0:
            return Trunc(self.truncation_degree, [a * 0 for a in self.coeffs])
        if not isinstance(other, Trunc):
            return self * Trunc.constant(self.truncation_degree, other)
        assert other.truncation_degree == self.truncation_degree
        zero = self.coeffs[0] * other.coeffs[0] * 0
        coeffs = []
        for k in range(self.truncation_degree + 1):
            total = zero
            for j in range(k + 1):
                total = total + self[j] * other[k - j]
            coeffs.append(total)
        return Trunc(self.truncation_degree, coeffs)

    @staticmethod
    def bound_high_order_product(self, other):
        assert isinstance(other, Trunc)
        assert other.truncation_degree == self.truncation_degree
        zero = self.coeffs[0] * other.coeffs[0] * 0
        total = zero
        for d in range(self.truncation_degree + 1, self.truncation_degree * 2 + 1):
            coeff = zero
            for k in range(min(d, self.truncation_degree) + 1):
                if 0 <= d - k <= self.truncation_degree:
                    coeff = coeff + self[k] * other[d - k]
            total = total + abs(coeff)
        return total

    @staticmethod
    def constant(n, a):
        zero = a * 0
        coeffs = [a] + [zero] * n
        return Trunc(n, coeffs)

    @staticmethod
    def basis_element(n, k, one):
        zero = one * 0
        coeffs = [one if j == k else zero for j in range(n + 1)]
        return Trunc(n, coeffs)

    def __rmul__(self, other):
        return self * other

    def __pow__(self, p):
        assert isinstance(p, int)
        assert p >= 0
        zero = self.coeffs[0] * 0
        one = self.coeffs[0] ** 0
        # the following relies on immutability of coeffs?
        # might be advisable to ensure separate zeros created.
        coeffs = [one] + [zero] * self.truncation_degree
        prod = Trunc(self.truncation_degree, coeffs)
        power = self
        while p > 0:
            if (p % 2) == 1:
                prod = prod * power
            power = power * power
            p = p // 2
        return prod

    def __call__(self, x):
        # horner algorithm is much more efficient
        #return sumup((x ** k)*a for k, a in enumerate(self.coeffs))
        ans = None
        for a in reversed(self.coeffs):
            if ans is None:
                ans = a
            else:
                ans = x*ans + a
        return ans

    def __abs__(self):
        return sumup(map(abs, self.coeffs))

    def abs_star(self):
        return sumup(map(abs, self.coeffs[1:]))

    def diff(self):
        if len(self) == 1:
            return Trunc(self.truncation_degree, [self[0] * 0])
        else:
            return Trunc(self.truncation_degree, [k * a for k, a in islice(enumerate(self), 1, None)])

