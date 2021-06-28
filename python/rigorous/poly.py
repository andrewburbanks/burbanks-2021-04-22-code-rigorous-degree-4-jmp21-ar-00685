def sumup(seq):
    elts = iter(seq)
    total = next(elts)
    for x in elts:
        total = total + x
    return total


from itertools import zip_longest, islice


def trim(coeffs):
    # this should be fairly fast
    # this works for intervals by using 0*coeffs[k] below
    for k in reversed(range(len(coeffs))):
        if not coeffs[k] == 0 * coeffs[k]:
            break
    return coeffs[:max(1, k + 1)]


class Poly:
    """
    We represent the zero polynomial by the degree zero Poly1d([0]).
    This ensures that at least one coefficient is present.
    This coefficient is used in various places to produce the relevant unit.

    The safest way to test if a polynomial is zero is whether self*0 == self.
    """

    def __init__(self, coeffs):
        assert len(coeffs) >= 1
        self.coeffs = list(trim(coeffs))

    @property
    def degree(self):
        return len(self.coeffs) - 1

    def __repr__(self):
        return 'Poly(%r)' % (self.coeffs)

    def __len__(self):
        return len(self.coeffs)

    def __iter__(self):
        yield from self.coeffs

    def __getitem__(self, k):
        # __iter__ is used for iteration, if present, otherwise
        # we fall back to __getitem__.
        # we must restrict k otherwise there is no stopiteration
        # the len is being ignored during iteration unless we
        # override iter directly.
        if k >= len(self.coeffs) or k < 0:
            raise IndexError
        return self.coeffs[k]

    def __eq__(self, other):
        # trust that other types know how to compare with int zero
        return all(a == b for a, b in zip_longest(self.coeffs, other.coeffs, fillvalue=0))

    def __pos__(self):
        return Poly([+a for a in self.coeffs])

    def __neg__(self):
        return Poly([-a for a in self.coeffs])

    def __add__(self, other):
        assert isinstance(other, Poly)
        return Poly([a + b for a, b in zip_longest(self.coeffs, other.coeffs, fillvalue=0)])

    def __sub__(self, other):
        assert isinstance(other, Poly)
        return Poly([a - b for a, b in zip_longest(self.coeffs, other.coeffs, fillvalue=0)])

    def __mul__(self, other):
        if isinstance(other, int) and other == 0:
            zero = self.coeffs[0] * 0
            return Poly([zero])
        if not isinstance(other, Poly):
            return self * Poly.constant(other)
        zero = self.coeffs[0] * other.coeffs[0] * 0
        coeffs = []
        for k in range(self.degree + other.degree + 1):
            total = zero
            for j in range(k + 1):
                if j < len(self) and 0 <= k - j < len(other):
                    total = total + self[j] * other[k - j]
            coeffs.append(total)
        return Poly(coeffs)

    @staticmethod
    def truncated_product(self, other, n):
        if isinstance(other, int) and other == 0:
            zero = self.coeffs[0] * 0
            return Poly([zero])
        if not isinstance(other, Poly):
            return self * Poly.constant(other)
        zero = self.coeffs[0] * other.coeffs[0] * 0
        coeffs = []
        for k in range(min(n, self.degree + other.degree) + 1):
            total = zero
            for j in range(k + 1):
                if j < len(self) and 0 <= k - j < len(other):
                    total = total + self[j] * other[k - j]
            coeffs.append(total)
        return Poly(coeffs)

    @staticmethod
    def constant(a):
        return Poly([a])

    @staticmethod
    def basis_element(n, k, one):
        zero = one * 0
        coeffs = [one if j == k else zero for j in range(n + 1)]
        return Poly(coeffs)

    def __rmul__(self, other):
        return self * other

    def __pow__(self, p):
        assert isinstance(p, int)
        assert p >= 0
        zero = self.coeffs[0] * 0
        one = self.coeffs[0] ** 0
        # the following relies on immutability of coeffs?
        # might be advisable to ensure separate zeros created.
        prod = Poly([one])
        power = self
        while p > 0:
            if (p % 2) == 1:
                prod = prod * power
            power = power * power
            p = p // 2
        return prod

    def __call__(self, x):
        # zero = x*0
        # total = zero
        # for k, a in enumerate(self.coeffs):
        #    total = total + a*x**k
        # return total
        return sumup(a * x ** k for k, a in enumerate(self.coeffs))

    def __abs__(self):
        return sumup(map(abs, self.coeffs))

    def truncate(self, truncation_degree):
        degree = min(truncation_degree, self.degree)
        return Poly(self.coeffs[:degree + 1])

    def split(self, truncation_degree):
        zero = self.coeffs[0] * 0
        P_lo = self.truncate(truncation_degree)
        P_hi = self - P_lo
        return P_lo, P_hi

    def diff(self):
        if len(self) == 1:
            return Poly([self[0] * 0])
        else:
            return Poly([k * a for k, a in islice(enumerate(self), 1, None)])

    @staticmethod
    def read(s):
        line = next(s)
        assert line.strip() == 'Poly1', line
        line = next(s)
        n = int(line.strip())
        coeffs = []
        for k in range(0, n + 1):
            line = next(s)
            sj, sa = line.strip().split()
            j, a = int(sj), float(sa)
            assert j == k
            coeffs.append(a)
        line = next(s)
        assert line.strip() == 'end.', line
        return Poly(coeffs)
