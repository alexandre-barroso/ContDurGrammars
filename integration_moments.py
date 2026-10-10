"""Three grid-centred moment arrays for compositional integral enclosures.

Signed convolutions combine disjoint children. Nonlinear factor updates retain
mass and first/second moments under derivative and support bounds. See METHOD.md.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F

from flint import arb, arb_poly
from integration_bounds import Potential, ball, tail_envelope


def interval(lo, hi):
    if lo > hi:
        raise ArithmeticError('reversed enclosing endpoints')
    return lo.union(hi)


def nonnegative_intersection(x, upper):
    return x.intersection(interval(arb(0), upper))


def leaf_moments(potential, lo, hi, center):
    """Exact Gaussian primitive formulas, evaluated using enclosing balls."""
    lo, hi, center = F(lo), F(hi), F(center)
    cuts = [lo, *(x for x in potential.knots if lo < x < hi), hi]
    totals = [arb(0), arb(0), arb(0)]
    for left, right in zip(cuts, cuts[1:]):
        a,b,c = potential.coefficients((left+right)/2)
        if a == 0:
            if b != 0:
                raise ValueError('unsupported linear-only leaf branch')
            scale = (-ball(c)).exp()
            moments = [scale*ball((right-center)**(j+1)-(left-center)**(j+1))/(j+1)
                       for j in range(3)]
        else:
            m0 = potential.gaussian_piece(left,right,(a,b,c))
            mu = b/a
            gl = (-ball(a*left*left-2*b*left+c)).exp()
            gr = (-ball(a*right*right-2*b*right+c)).exp()
            centered1 = (gl-gr)/ball(2*a)
            centered2 = (m0-ball(right-mu)*gr+ball(left-mu)*gl)/ball(2*a)
            shift = ball(mu-center)
            moments = [m0, centered1+shift*m0,
                       centered2+2*shift*centered1+shift*shift*m0]
        totals = [x+y for x,y in zip(totals,moments)]
    radius = ball(max(abs(lo-center),abs(hi-center)))
    # Bounds follow from positive measure and its known support, not a repair
    # of an inconsistent calculation: every intersection must be nonempty.
    totals[0] = nonnegative_intersection(totals[0],totals[0].upper().max(0))
    totals[1] = totals[1].intersection(arb(0,(radius*totals[0].upper()).upper()))
    totals[2] = nonnegative_intersection(totals[2],(radius**2*totals[0].upper()).upper())
    return totals


def factor_derivative_bounds(potential, lo, hi, center):
    """Bounds for |g'| and the Lipschitz constant of g' across C1 knots."""
    lo, hi, center = F(lo), F(hi), F(center)
    cuts = [lo, *(x for x in potential.knots if lo < x < hi), hi]
    k1 = k2 = arb(0)
    for left,right in zip(cuts,cuts[1:]):
        a,b,c = potential.coefficients((left+right)/2)
        minimizer = min(right,max(left,b/a)) if a else left
        max_factor = (-ball(a*minimizer**2-2*b*minimizer+c)).exp().upper()
        dl, dr = 2*(a*left-b), 2*(a*right-b)
        derivative_max = max(abs(dl),abs(dr))
        square_min = F(0) if dl <= 0 <= dr else min(dl*dl,dr*dr)
        square_max = max(dl*dl,dr*dr)
        second_max = max(abs(square_min-2*a),abs(square_max-2*a))
        k1 = k1.max((ball(derivative_max)*max_factor).upper())
        k2 = k2.max((ball(second_max)*max_factor).upper())
    a,b,_ = potential.coefficients(center)
    g0 = (-potential.value(center)).exp()
    g1 = -ball(2*(a*center-b))*g0
    return g0,g1,k1.upper(),k2.upper()


@dataclass
class MomentMessage:
    moments: list[list]
    leaves: int
    origin: F


def combine(left, right):
    a = [arb_poly(x) for x in left.moments]
    b = [arb_poly(x) for x in right.moments]
    polynomials = [a[0]*b[0], a[1]*b[0]+a[0]*b[1],
                   a[2]*b[0]+2*a[1]*b[1]+a[0]*b[2]]
    length = len(left.moments[0])+len(right.moments[0])-1
    return MomentMessage([[p[k] for k in range(length)] for p in polynomials],
                         left.leaves+right.leaves,left.origin+right.origin)


def multiply_factor(moments, potential, center, radius):
    center, radius = F(center), F(radius)
    lo,hi = center-radius,center+radius
    m0,m1,m2 = moments
    r = ball(radius)
    # Convolution inputs enclose true nonnegative 0th/2nd moments. The balls
    # themselves may extend below zero, so use an exact nonnegative upper cap.
    m0 = nonnegative_intersection(m0,m0.upper().max(0))
    cap0 = m0.upper()
    m1 = m1.intersection(arb(0,(r*cap0).upper()))
    m2 = nonnegative_intersection(m2,(r*r*cap0).upper())
    cap2 = m2.upper().max(0)
    g0,g1,k1,k2 = factor_derivative_bounds(potential,lo,hi,center)
    a0 = g0*m0+g1*m1+arb(0,(k2*cap2/2).upper())
    a1 = g0*m1+g1*m2+arb(0,(k2*r*cap2/2).upper())
    a2 = g0*m2+arb(0,(k1*r*cap2).upper())
    lower,upper = potential.factor_bounds(lo,hi)
    factor = interval(lower,upper)
    a0 = a0.intersection(factor*m0)
    a0 = nonnegative_intersection(a0,a0.upper().max(0))
    a1 = a1.intersection(arb(0,(r*a0.upper()).upper()))
    a2 = a2.intersection(factor*m2)
    a2 = nonnegative_intersection(a2,(r*r*a0.upper()).upper())
    if not all(x.is_finite() for x in (a0,a1,a2)):
        raise ArithmeticError('nonfinite moment enclosure')
    return a0,a1,a2


def enclose_moments(node,h,lo,hi,cache):
    h, lo, hi = F(h), F(lo), F(hi)
    if not node.children:
        key = (node.potential,h,lo,hi)
        if key not in cache:
            cells = (hi-lo)/h
            if cells.denominator != 1:
                raise ValueError('mesh and bounds must align')
            values = [leaf_moments(node.potential,lo+k*h,lo+(k+1)*h,lo+(F(k)+F(1,2))*h)
                      for k in range(int(cells))]
            cache[key] = MomentMessage([list(x) for x in zip(*values)],1,lo)
        return cache[key]
    children = [enclose_moments(c,h,lo,hi,cache) for c in node.children]
    while len(children)>1:
        children = [combine(children[i],children[i+1]) if i+1<len(children) else children[i]
                    for i in range(0,len(children),2)]
    message = children[0]
    radius = message.leaves*h/2
    result = [[],[],[]]
    for k,moments in enumerate(zip(*message.moments)):
        center = message.origin+k*h+radius
        values = multiply_factor(moments,node.potential,center,radius)
        for row,value in zip(result,values):
            row.append(value)
    return MomentMessage(result,message.leaves,message.origin)


def partition_enclosure_moments(node,h,hi):
    h, hi = F(h), F(hi)
    message = enclose_moments(node,h,F(0),hi,{})
    mass = sum(message.moments[0],arb(0))
    lower = mass.lower().max(0)
    tail = tail_envelope(node,hi)
    upper = (mass.upper()+tail).upper()
    width = (upper.log()-lower.log()).upper() if lower > 0 else arb.pos_inf()
    return lower,upper,tail,width,len(message.moments[0])
