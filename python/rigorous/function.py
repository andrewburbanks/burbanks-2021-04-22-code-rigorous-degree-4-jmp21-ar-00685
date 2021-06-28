from .interval import Interval, RR, str_prec
from .rectangle import Rectangle
from .ball import Ball
from .trunc import Trunc

ArgTypes = (Rectangle, Interval, RR, float, int)

class Domain:

    """
    A disc domain.  The centre and radius can be non-rigorous or rigorous.
    """

    def __init__(self, centre, radius):
        self.c = centre
        self.r = radius

    def to_unit(self, x):
        u = (x - self.c) * (1 / self.r)
        return u

    def from_unit(self, u):
        x = self.c + u * self.r
        return x

    def __eq__(self, other):
        assert isinstance(other, Domain)
        return self.c == other.c and self.r == other.r

    def __repr__(self):
        return f'Domain(centre={self.c}, radius={self.r})'


class Function:

    """
    A Function object represents a power series expanded on some domain.
    """

    def __init__(self, dom, func):
        self.dom = dom
        self.func = func

    def __eq__(self, other):
        assert isinstance(other, Function)
        return self.dom == other.dom and self.func == other.func

    # some ugly code here, testing the types. it should be possible to unify
    # these by placing some assumptions on the arguments, i.e., that they
    # support the relevant operations.

    def promote(self, x):
        """
        When promoting any objects, we assume the default domain (0, 1).
        It is important to do this, rather than assuming self.dom, in
        order that we catch errors where we've neglected to take domain
        into account.
        """
        
        if isinstance(x, Function):
            return x
        if isinstance(x, Ball):
            return Function(Domain(0.0, 1.0), x)
        if isinstance(x, Trunc):
            return Function(Domain(0.0, 1.0), Ball(x, 0.0, 0.0))
        if isinstance(x, ArgTypes):
            return x
        raise NotImplementedError

    def __call__(self, x):
        x = self.promote(x)
        u = self.dom.to_unit(x)
        if isinstance(x, Function):
            return Function(x.dom, self.func(u.func))
        if isinstance(x, ArgTypes):
            return self.func(u)
        raise NotImplementedError

    def __mul__(self, x):
        x = self.promote(x)
        if isinstance(x, Function):
            assert self.dom == x.dom, 'Domain mismatch!'
            return Function(self.dom, self.func * x.func)
        if isinstance(x, ArgTypes):
            return Function(self.dom, self.func * x)
        raise NotImplementedError

    def __truediv__(self, x):
        if isinstance(x, float) or isinstance(x, int):
            x = Interval(x, x)
        assert isinstance(x, Interval)
        return self*(1/x)
    
    def __pow__(self, p):
        assert isinstance(p, int)
        assert p >= 0
        one = Function(self.dom, self.func ** 0)
        prod = one
        power = self
        while p > 0:
            if (p % 2) == 1:
                prod = prod * power
            power = power * power
            p = p // 2
        return prod

    def __rmul__(self, x):
        return self * x

    def __add__(self, x):
        x = self.promote(x)
        if isinstance(x, Function):
            assert self.dom == x.dom, 'Domain mismatch!'
            return Function(self.dom, self.func + x.func)
        if isinstance(x, ArgTypes):
            return Function(self.dom, self.func + x)
        raise NotImplementedError

    def __radd__(self, x):
        return self + x

    def __sub__(self, x):
        x = self.promote(x)
        if isinstance(x, Function):
            assert self.dom == x.dom, 'Domain mismatch!'
            return Function(self.dom, self.func - x.func)
        if isinstance(x, ArgTypes):
            return Function(self.dom, self.func - x)
        raise NotImplementedError

    def __rsub__(self, x):
        return -self+x

    def __pos__(self):
        return self

    def __neg__(self):
        return Function(self.dom, -self.func)

    def __abs__(self):
        return abs(self.func)

    def __repr__(self):
        return f'Function(domain={self.dom}, func={self.func})'

    def identity(self):
        i1 = Interval(1, 1)
        if isinstance(self.func, Ball):
            N = self.func.truncation_degree
            x = Function(self.dom, Ball(Trunc(N, [self.dom.c*i1, self.dom.r*i1]), 0, 0))
        elif isinstance(self.func, Trunc):
            N = self.func.truncation_degree
            x = Function(self.dom, Trunc(N, [self.dom.c*i1, self.dom.r*i1]), 0, 0)
        else:
            raise NotImplementedError
        return x

    def zero(self):
        return self.identity()*0

    # we need a dcomp member function for derivative compose.
    def diff_compose(self, other):
        other = self.promote(other)
        if isinstance(self.func, Trunc):
            u = self.dom.to_unit(other)
            return Function(other.dom, self.func.diff()(u.func)*(1 / self.dom.r)) #u.func?
        assert isinstance(self.func, Ball)
        u = self.dom.to_unit(other)
        return Function(other.dom, self.func.diff_compose(u.func) * (1 / self.dom.r)) #u.func?

    def __str__(self):
        domain, ball = self.dom, self.func
        s0 = f'Function Ball:'
        s0 += f' c={domain.c}, r={domain.r}, N={ball.truncation_degree} ({str_prec} sf shown)'
        s1 = '\n'.join(f'{k:0>3}: '+str(a) for k, a in enumerate(ball.P))
        s2 = f'  H: {str(ball.H)}'
        s3 = f'  G: {str(ball.G)}'
        return '\n'.join([s0, s1, s2, s3])

    def basis_element(self, k):
        assert 0 <= k <= self.func.truncation_degree+1
        e_k = Ball.basis_element(self.func.truncation_degree, k, Interval(1))
        return Function(self.dom, e_k)

    def high_order_ball(self, rho):
        N = self.func.truncation_degree
        x = Function(self.dom, Ball(Trunc(N, [Interval(0)]), rho, 0))
        return x

    def general_ball(self, rho):
        N = self.func.truncation_degree
        x = Function(self.dom, Ball(Trunc(N, [Interval(0)]), 0, rho))
        return x
