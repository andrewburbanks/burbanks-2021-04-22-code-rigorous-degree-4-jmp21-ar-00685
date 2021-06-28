from .interval import Interval, RR

class Rectangle:

    def __init__(self, re, im=None):
        if isinstance(re, (RR, float, int)):
            re = Interval(re, re)
        if im is None:
            im = Interval(0, 0)
        if isinstance(im, (RR, float, int)):
            im = Interval(im, im)
        assert isinstance(re, Interval)
        assert isinstance(im, Interval)
        self.re = re
        self.im = im

    def __repr__(self):
        return f'Rectangle(re={self.re}, im={self.im})'

    # inspectors
    
    def real(self):
        return self.re

    def imag(self):
        return self.im

    # unary arithmetic operators

    def __pos__(self):
        re = self.re
        im = self.im
        return Rectangle(re, im)

    def __neg__(self):
        re = -self.re
        im = -self.im
        return Rectangle(re, im)

    # binary arithmetic operators

    def __add__(self, other):
        if isinstance(other, (Interval, RR, float, int)):
            other = Rectangle(other)
        assert isinstance(other, Rectangle)
        re = self.re + other.re
        im = self.im + other.im
        return Rectangle(re, im)

    def __sub__(self, other):
        if isinstance(other, (Interval, RR, float, int)):
            other = Rectangle(other)
        assert isinstance(other, Rectangle)
        re = self.re - other.re
        im = self.im - other.im
        return Rectangle(re, im)

    def __mul__(self, other):
        if isinstance(other, (Interval, RR, float, int)):
            other = Rectangle(other)
        assert isinstance(other, Rectangle)
        re = self.re * other.re - self.im * other.im
        im = self.re * other.im + self.im * other.re
        return Rectangle(re, im)

    def __truediv__(self, other):
        if isinstance(other, (Interval, RR, float, int)):
            return Rectangle(self.re/other, self.im/other)
        assert isinstance(other, Rectangle)
        assert not (0 in other.re and 0 in other.im)
        denom = other.re**2 + other.im**2
        re = (self.re*other.re + self.im*other.im)/denom
        im = (self.im*other.re - self.re*other.im)/denom
        return Rectangle(re, im)

    def __pow__(self, p):
        if isinstance(p, (RR, float, int)) and p == 0:
            return Rectangle(1)
        assert isinstance(p, int)
        prod = Rectangle(1)
        ppow = self
        while p > 0:
            if p%2 == 1:
                prod = prod * ppow
            p = p // 2
            ppow = ppow * ppow
        return prod

    # right- versions of binary arithmetic operators

    def __radd__(self, other):
        return self+other

    def __rsub__(self, other):
        return -(self-other)
    
    def __rmul__(self, other):
        return self*other

    def __rtruediv__(self, other):
        raise NotImplementedError

    # modulus

    def __abs__(self):
        return (self.re**2 + self.im**2)**0.5
