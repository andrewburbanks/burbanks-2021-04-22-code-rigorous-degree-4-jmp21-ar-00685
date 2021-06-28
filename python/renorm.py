import sys

from rigorous.poly import Poly
from rigorous.trunc import Trunc
from rigorous.interval import Interval, RR
from rigorous.rectangle import Rectangle
from rigorous.covering import cover_circle
from rigorous.ball import Ball
from rigorous.function import Domain, Function

# we work in the symmetry-reduced subspace
# unlike lanford, we keep the constant term in our series

d = 2 #override as needed

def Q(G):
    return G**d

def Q_dash(G):
    return G**(d-1)*d

def T_even(G):
    assert isinstance(G, Function)

    x = G.identity()
    alpha = 1/G(1)
    return G(Q(G(x/alpha**d)))*alpha

def DT_even(G, H):
    assert isinstance(G, Function)
    assert isinstance(H, Function)

    x = G.identity()
    a = G(1)
    alpha = 1/a
    dG = G.diff_compose
    deltaT1 = G(Q(G(x*a**d)))*(-alpha**2*H(1))
    deltaT2 = H(Q(G(x*a**d)))*alpha
    deltaT3 = dG(Q(G(x*a**d)))*alpha*Q_dash(G(x*a**d))*H(x*a**d)
    deltaT4 = dG(Q(G(x*a**d)))*alpha*Q_dash(G(x*a**d))*dG(x*a**d)*Q_dash(a)*x*H(1)
    deltaT = deltaT1 + deltaT2 + deltaT3 + deltaT4
    return deltaT

def DT(G):
    assert isinstance(G, Function)
    
    # compute constants
    a = G(1)
    alpha = 1/a

    # compute chain in TG
    X = G.identity()
    Xad = X*a**d
    GXad = G(Xad)
    QGXad = Q(GXad)
    GQGXad = G(QGXad)
    
    # compute derivatives
    Gdash = G.diff_compose
    GdashXad = Gdash(Xad)
    QdashGXad = Q_dash(GXad)
    GdashQGXad = Gdash(QGXad)

    # compute all prefactors for frechet derivative
    A = GQGXad*(-alpha**2)
    B = alpha
    C = GdashQGXad*QdashGXad*alpha
    D = C*Q_dash(a)*X
    
    def DTG(H):
        assert isinstance(H, Function)
        H1 = H(1)
        return A*H1 + H(QGXad)*B + C*H(Xad) + D*H1
    
    return DTG

#def g_even_dx_at_fixed_point(G):
#    assert isinstance(G, Function)
#
#    X = G.identity()
#    a = G(1)
#    dG = G.diff_compose
#    a2X = X*a**2
#    Ga2X = G(a2X)
#    QGa2X = Q(Ga2X)
#    return dG(QGa2X)*2*Ga2X*dG(a2X)*a

# load domain and coefficients of a truncated power series
def import_poly_with_domain(filename):
    with open(filename, 'r') as lines:
        c = RR(next(lines))
        r = RR(next(lines))
        d = int(next(lines))
        coeffs = []
        for k in range(d+1):
            coeffs.append(RR(next(lines)))
    return c, r, coeffs

# load matrices
def import_matrix(filename):
    with open(filename, 'r') as f:
        m = int(next(f))
        n = int(next(f))
        assert m == n
        rows = []
        for k in range(m):
            row = list(map(lambda x: Interval(RR(x), RR(x)), next(f).strip().split(' ')))
            assert len(row) == n
            rows.append(row)
    return rows

# we had a matrix Delta that approximates the derivative of T
# we inverted (Delta - I) to produce a matrix Lambda for the newton method.
# we now define a linear operator L in our newton method
# the action of this operator
# on the polynomial part of the space is Lambda=(Delta-I)^{-1}
# on the high-order part of the space is -I
def apply_lambda_operator(Lambda, F):
    trunc_degree = F.func.truncation_degree
    assert isinstance(Lambda, list)
    assert len(Lambda) == len(Lambda[0]) == trunc_degree+1
    assert isinstance(F, Function)
    coeffs = []
    for row in Lambda:
        coeff = sum(a*b for a, b in zip(row, F.func.P))
        # worst-case scenario, G is in a single polynomial coefficient
        coeff = coeff + Interval.hull(row) * Interval(-F.func.G.hi, F.func.G.hi)
        coeffs.append(coeff)
    P = Trunc(trunc_degree, coeffs)
    # action on high-order part of the space is trivial
    H = F.func.H + F.func.G # because we have -I in the H/H quadrant
    # action on general part of space, bounded by worst case
    G = 0 #F.func.G # NOT needed in addition to the above; overkill
    # note -f_H bounded in size by +H.
    return Function(F.dom, Ball(P, H, G))

