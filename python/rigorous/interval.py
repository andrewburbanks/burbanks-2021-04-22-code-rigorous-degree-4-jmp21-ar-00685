from decimal import Decimal, getcontext
from decimal import localcontext, ROUND_CEILING, ROUND_FLOOR, Inexact

RR = Decimal

str_prec = 15

def set_precision(prec):
    getcontext().prec = prec
    
def set_output_precision(prec):
    global str_prec
    str_prec = prec
    
def quantize(x):
    """
    
    Safely convert interval to a simpler one in which only the final dp differs.
    This is used to indicate which digits of the number are known for certain.
    
    """
    rho = RR(10)**(int(abs(x.hi-x.lo).log10()))
    a = x.lo.quantize(rho, rounding=ROUND_FLOOR)
    b = x.hi.quantize(rho, rounding=ROUND_CEILING)
    return Interval(a, b)

def add_up(a, b):
    with localcontext() as ctx:
        ctx.rounding = ROUND_CEILING
        return a+b

def add_dn(a, b):
    with localcontext() as ctx:
        ctx.rounding = ROUND_FLOOR
        return a+b

def mul_up(a, b):
    with localcontext() as ctx:
        ctx.rounding = ROUND_CEILING
        return a*b

def mul_dn(a, b):
    with localcontext() as ctx:
        ctx.rounding = ROUND_FLOOR
        return a*b

def sqrt_up(a):
    if a == 0:
        return RR(0)
    with localcontext() as ctx:
        ctx.rounding = ROUND_CEILING
        hi = a.sqrt().next_plus()
        # ROUND_NEAR overrides selected mode for sqrt, log, etc.
        ctx.rounding = ROUND_FLOOR
        assert a <= hi*hi
        return hi

def sqrt_dn(a):
    if a == 0:
        return RR(0)
    with localcontext() as ctx:
        ctx.rounding = ROUND_FLOOR
        lo = a.sqrt().next_minus()
        # ROUND_NEAR overrides selected mode for sqrt, log, etc.
        ctx.rounding = ROUND_CEILING
        assert lo*lo <= a
        return lo

def div_up(a, b):
    with localcontext() as ctx:
        ctx.rounding = ROUND_CEILING
        return a/b

def div_dn(a, b):
    with localcontext() as ctx:
        ctx.rounding = ROUND_FLOOR
        return a/b

def square(x):
    assert isinstance(x, Interval)
    if x.lo <= 0 <= x.hi:
        lo = 0
    else:
        lo = max(0, min(mul_dn(x.lo, x.lo), mul_dn(x.hi, x.hi)))
    hi = max(mul_up(x.lo, x.lo), mul_up(x.hi, x.hi))
    return Interval(lo, hi)

def int_power(x, p):
    assert isinstance(p, int)
    if p < 0:
        return Interval(1, 1)/power(x, -p)
    if p == 0:
        return Interval(1, 1)
    elif p == 1:
        return Interval(x.lo, x.hi)
    xpow = x
    prod = None #avoid multiplication by 1
    while p > 0:
        if p%2 == 1:
            if prod is None:
                prod = xpow
            else:
                prod = prod * xpow
        xpow = square(xpow)
        p = p//2
    return prod

class Interval:

    def __init__(self, lo, hi=None):
        if hi is None:
            if isinstance(lo, Interval):
                lo, hi = lo.lo, lo.hi
            else:
                hi = lo
        self.lo = RR(lo)
        self.hi = RR(hi)
        assert self.lo <= self.hi

    def __repr__(self):
        return 'Interval(%r, %r)' % (self.lo, self.hi)

    def __eq__(self, other):
        assert isinstance(other, Interval)
        return self.lo == other.lo and self.hi == other.hi

    def __neg__(self):
        lo = -self.hi
        hi = -self.lo
        return Interval(lo, hi)

    def __pos__(self):
        lo = +self.lo
        hi = +self.hi
        return Interval(lo, hi)

    def __add__(self, other):
        if isinstance(other, (int, float, RR)):
            other = Interval(other, other)

        assert isinstance(other, Interval)
        with localcontext() as ctx:
            ctx.rounding = ROUND_FLOOR
            lo = self.lo + other.lo
            ctx.rounding = ROUND_CEILING
            hi = self.hi + other.hi
        return Interval(lo, hi)

    def __sub__(self, other):
        """Ensure that p-q == p+(-q) by defining it that way."""
        if isinstance(other, (int, float, RR)):
            other = Interval(other, other)

        assert isinstance(other, Interval)
        return self + (-other)

    def __mul__(self, other):
        if isinstance(other, int) and other == 0:
            return Interval(self.lo * 0, self.hi * 0)
        if isinstance(other, (int, float, RR)):
            other = Interval(other, other)

        assert isinstance(other, Interval)
        with localcontext() as ctx:
            ctx.rounding = ROUND_FLOOR
            lo = min(self.lo * other.lo,
                     self.lo * other.hi,
                     self.hi * other.lo,
                     self.hi * other.hi)
            ctx.rounding = ROUND_CEILING
            hi = max(self.lo * other.lo,
                     self.lo * other.hi,
                     self.hi * other.lo,
                     self.hi * other.hi)
        return Interval(lo, hi)

    def inverse(self):
        if 0 in self:
            raise ZeroDivisionError
        lo = div_dn(1, self.hi)
        hi = div_up(1, self.lo)
        return Interval(lo, hi)

    def __truediv__(self, other):
        if isinstance(other, (int, float, RR)):
            other = Interval(other, other)
        assert isinstance(other, Interval)
        return self * other.inverse()

    def __radd__(self, other):
        return self+other

    def __rsub__(self, other):
        return -(self-other)

    def __rmul__(self, other):
        return self*other

    def __rtruediv__(self, other):
        return self.inverse()*other

    def __pow__(self, other):
        assert isinstance(other, (int, float, RR))
        if other < 0:
            return Interval(1)/self**(-other)
        if other == 0:
            return Interval(1, 1)
        if other == 1:
            return Interval(self.lo, self.hi)
        if not isinstance(other, int):
            return (self.log()*RR(other)).exp()
        return int_power(self, other)

    def __abs__(self):
        a, b = abs(self.lo), abs(self.hi)
        if 0 in self:
            return Interval(a * 0, max(a, b))
        else:
            return Interval(min(a, b), max(a, b))

    def __lt__(self, other):
        if isinstance(other, (int, float, RR)):
            other = Interval(other, other)
        assert isinstance(other, Interval)
        return self.hi < other.lo

    def __le__(self, other):
        if isinstance(other, (int, float, RR)):
            other = Interval(other, other)
        assert isinstance(other, Interval)
        return self.hi <= other.lo

    def __ge__(self, other):
        if isinstance(other, (int, float, RR)):
            other = Interval(other, other)
        assert isinstance(other, Interval)
        return self.lo >= other.hi

    def exp(self):
        # the ibm decimal standard promises at most 1ulp error
        lo = self.lo.exp().next_minus()
        hi = self.hi.exp().next_plus()
        return Interval(lo, hi)

    def log(self):
        # the ibm decimal standard promises at most 1ulp error
        assert 0 <= self.lo
        lo = self.lo.ln().next_minus()
        hi = self.hi.ln().next_plus()
        return Interval(lo, hi)

    def __contains__(self, other):
        if isinstance(other, Interval):
            return self.lo <= other.lo and other.hi <= self.hi
        return self.lo <= other <= self.hi

    @staticmethod
    def hull(intervals):
        lo = min(a.lo for a in intervals)
        hi = max(a.hi for a in intervals)
        return Interval(lo, hi)

    def __str__(self):
        with localcontext() as ctx:
            ctx.prec = str_prec
            if self.lo == 0:
                lo = 0
            else:
                ctx.rounding = ROUND_FLOOR
                lo = RR(self.lo+0)
            if self.hi == 0:
                hi = 0
            else:
                ctx.rounding = ROUND_CEILING
                hi = RR(self.hi+0)
        assert lo <= self.lo <= self.hi <= hi
        return f'[{lo:+1.{str_prec-1}E}, {hi:+1.{str_prec-1}E}]'
    
def Pub(x):
    if isinstance(x, (int, float, RR)):
        return Interval(0, abs(x))
    assert isinstance(x, Interval)
    return Interval(0, abs(x).hi)
