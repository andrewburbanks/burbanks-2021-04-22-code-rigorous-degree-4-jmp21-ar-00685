from rigorous.interval import RR, Interval
from rigorous.trunc import Trunc
from rigorous.ball import Ball
from rigorous.function import Function

class RMatrix:

    """
    Matrices of intervals.
    """

    @staticmethod
    def promote(x):
        if isinstance(x, (int, float)):
            return Interval(RR(x))
        if isinstance(x, RR):
            return Interval(x)
        if isinstance(x, Interval):
            return x
        raise NotImplementedError

    @staticmethod
    def zeros(size):
        m, n = size
        a = [[Interval(0) for _ in range(n)] for _ in range(m)]
        return RMatrix(a)
    
    @staticmethod
    def eye(n):
        A = RMatrix.zeros((n, n))
        for k in range(n):
            A.elts[k][k] = Interval(1)
        return A
    
    @staticmethod
    def basis_element(n, k):
        a = [[Interval(0)] for j in range(n)]
        a[k][0] = Interval(1)
        return RMatrix(a)
        
    def __init__(self, a):
        self.elts = [[None for _ in a[0]] for _ in a]
        for j in range(self.m):
            for k in range(self.n):
                self.elts[j][k] = RMatrix.promote(a[j][k])
                
    def __repr__(self):
        return str(self)
    
    def __str__(self):
        return '['+',\n '.join(['['+', '.join(map(str, row))+']' for row in self.elts])+']'

    def copy(self):
        return RMatrix(self.elts)
    
    def row(self, j):
        return RMatrix([self.elts[j]])

    def col(self, k):
        return RMatrix([[row[k]] for row in self.elts])

    @property
    def m(self):
        return len(self.elts)
    @property
    def n(self):
        return len(self.elts[0])
    
    def shape(self):
        return (self.m, self.n)

    def T(self):
        return RMatrix(list(zip(*self.elts)))
    
    def column_norm(self, k):
        return sum([abs(self.elts[j][k]) for j in range(self.m)], Interval(0))
    
    def __abs__(self):
        return Interval.hull([self.column_norm(k) for k in range(self.n)])

    def __eq__(self, that):
        m, n = self.shape()
        o, p = that.shape()
        if not (m == o and n == p):
            return False
        return all(r == s for r, s in zip(self.elts, that.elts))
  
    @staticmethod
    def concatenate(A, B):
        m, n = A.shape()
        o, p = B.shape()
        assert m == o
        return RMatrix([r+s for r, s in zip(A.elts, B.elts)])

    @staticmethod
    def split(A, n):
        m, nn = A.shape()
        assert n < nn
        B = RMatrix([r[:n] for r in A.elts])
        C = RMatrix([r[n:] for r in A.elts])
        return B, C
    
    @staticmethod
    def gauss(A):
        """
        Gaussian elimation.
        Works for non-square matrices for use as solver and inverter.
        """
        A = A.copy()      
        m, n = A.shape()
        assert m <= n
        a = A.elts
        
        for k in range(m):
            # seek largest nonzero a[j][k] across j
            best, abs_best = k, abs(a[k][k])
            for i in range(k+1, m):
                if abs(a[i][k]) > abs_best:
                    best, abs_best = i, abs(a[i][k])

            # swap rows to position pivot
            a[k], a[best] = a[best], a[k]
            
            # divide row by pivot
            pivot = a[k][k]
            for j in range(k, n):
                a[k][j] /= pivot

            # clear zeros above and below
            for k_that in range(m):
                if k_that != k:
                    ratio = a[k_that][k]/a[k][k]
                    for j in range(k, n):
                        a[k_that][j] -= ratio*a[k][j]
        return A
   
    @staticmethod
    def det(A):
        A = A.copy()      
        m, n = A.shape()
        assert m <= n
        a = A.elts
        
        prod = Interval(1)

        for k in range(m):
            # seek largest nonzero a[j][k] across j
            best, abs_best = k, abs(a[k][k])
            for i in range(k+1, m):
                if abs(a[i][k]) > abs_best:
                    best, abs_best = i, abs(a[i][k])

            # swap rows to position pivot
            if k != best:
                a[k], a[best] = a[best], a[k]
                prod = -prod

            pivot = a[k][k]
            prod = prod * pivot

            # clear zeros below
            for k_that in range(k+1, m):
                if k_that != k:
                    ratio = a[k_that][k]/a[k][k]
                    for j in range(k, n):
                        a[k_that][j] -= ratio*a[k][j]
        return prod

    def inv(self):
        m, n = self.shape()
        assert m == n
        A = self
        I = RMatrix.eye(n)
        A_I = RMatrix.concatenate(A, I)
        I_Ainv = RMatrix.gauss(A_I)
        I_, Ainv = RMatrix.split(I_Ainv, n)
        return Ainv
    
    def __getitem__(self, index):
        j, k = index
        return self.elts[j][k]
    
    def __setitem__(self, index, x):
        j, k = index
        self.elts[j][k] = RMatrix.promote(x)
    
    def __mul__(self, that):
        if isinstance(that, (int, float, RR, Interval)):
            return RMatrix([[x*that for x in row] for row in self.elts])
        if isinstance(that, Function):
            return self(that)
        assert isinstance(that, RMatrix)
        m, n = self.shape()
        o, p = that.shape()
        assert n == o
        a, b = self.elts, that.elts
        c = [[sum((a[j][i]*b[i][k] for i in range(n)), Interval(0)) for k in range(p)] for j in range(m)]
        return RMatrix(c)
    
    def __rmul__(self, that):
        return self*that

    def __add__(self, that):
        if isinstance(that, (int, float, RR, Interval)):
            return RMatrix([[x+that for x in row] for row in self.elts])
        assert isinstance(that, RMatrix)
        m, n = self.shape()
        o, p = that.shape()
        assert m == o and n == p
        a, b = self.elts, that.elts
        c = [[a[j][k]+b[j][k] for k in range(p)] for j in range(m)]
        return RMatrix(c)
    
    def __radd__(self, that):
        return self+that
    
    def __pos__(self):
        return RMatrix([[+x for x in row] for row in self.elts])

    def __neg__(self):
        return RMatrix([[-x for x in row] for row in self.elts])

    def __sub__(self, that):
        return self+(-that)
    
    def __rsub__(self, that):
        return (-self)+that
    
    def __call__(self, F, high_order_coeff=1):
        """
        Apply the linear operator having this matrix as low-order
        part and +/-hI as the high-order part to a function ball.
        """
        assert isinstance(F, Function)
        trunc_degree = F.func.truncation_degree

        Lambda = self.elts
        assert isinstance(Lambda, list)
        assert len(Lambda) == len(Lambda[0]) == trunc_degree+1

        coeffs = []
        for row in Lambda:
            coeff = sum((a*b for a, b in zip(row, F.func.P)), Interval(0))

            # worst-case scenario, G is in a single polynomial coefficient
            coeff = coeff + Interval.hull(row) * Interval(-F.func.G.hi, F.func.G.hi)
            coeffs.append(coeff)
        P = Trunc(trunc_degree, coeffs)
        H = (F.func.H + F.func.G)*abs(high_order_coeff)
        G = 0
        return Function(F.dom, Ball(P, H, G))
