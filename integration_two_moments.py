"""Two-component normalizer calculation; see TWO_MOMENT.md.

Gaussian/Arb primitives and model constructors are shared with the three-component
algorithm. No second-moment array is computed or stored.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from flint import arb, arb_poly
from integration_bounds import ball, tail_envelope
from integration_moments import factor_derivative_bounds


def hull(lo, hi):
    if lo > hi:
        raise ArithmeticError("reversed interval")
    return lo.union(hi)


def mass_cap(value):
    return value.intersection(hull(arb(0), value.upper().max(0)))


def checked(values):
    if not all(x.is_finite() for x in values):
        raise ArithmeticError("nonfinite two-moment enclosure")
    return values


def leaf_two(potential, lo, hi, center):
    """Only Gaussian zeroth/first primitives; exact rational branching."""
    lo, hi, center = map(F, (lo, hi, center))
    if lo > hi:
        raise ValueError("reversed leaf interval")
    cuts = [lo, *(x for x in potential.knots if lo < x < hi), hi]
    zero, first = arb(0), arb(0)
    for left, right in zip(cuts, cuts[1:]):
        if left == right:
            continue
        a, b, c = potential.coefficients((left + right) / 2)
        if a == 0:
            if b != 0:
                raise ValueError("linear-only branch is outside the contract")
            scale = (-ball(c)).exp()
            zero += scale * ball(right-left)
            first += scale * ball(((right-center)**2-(left-center)**2) / 2)
        else:
            m0 = potential.gaussian_piece(left, right, (a, b, c))
            gl = (-ball(a*left*left-2*b*left+c)).exp()
            gr = (-ball(a*right*right-2*b*right+c)).exp()
            zero += m0
            first += (gl-gr)/ball(2*a) + ball(b/a-center)*m0
    radius = ball(max(abs(lo-center), abs(hi-center)))
    zero = mass_cap(zero)
    first = first.intersection(arb(0, (radius*zero.upper()).upper()))
    return checked((zero, first))


@dataclass
class TwoMessage:
    moments: list[list]
    leaves: int
    origin: F


def merge_two(left, right):
    a0, a1 = (arb_poly(row) for row in left.moments)
    b0, b1 = (arb_poly(row) for row in right.moments)
    values = (a0*b0, a1*b0+a0*b1)
    length = len(left.moments[0])+len(right.moments[0])-1
    return TwoMessage([[polynomial[k] for k in range(length)] for polynomial in values],
                      left.leaves+right.leaves, left.origin+right.origin)


def factor_two(moments, potential, center, radius):
    center, radius = F(center), F(radius)
    if radius < 0:
        raise ValueError("negative radius")
    m0, m1 = moments
    r = ball(radius)
    m0 = mass_cap(m0)
    upper0 = m0.upper()
    m1 = m1.intersection(arb(0, (r*upper0).upper()))
    g0, g1, k1, k2 = factor_derivative_bounds(
        potential, center-radius, center+radius, center)
    # A support-only bound, not a stored or computed second moment.
    support_second_upper = (r*r*upper0).upper()
    out0 = g0*m0+g1*m1+arb(0, (k2*support_second_upper/2).upper())
    out1 = g0*m1+arb(0, (k1*support_second_upper).upper())
    lower, upper = potential.factor_bounds(center-radius, center+radius)
    factor_range = hull(lower.max(0), upper.min(1))
    out0 = mass_cap(out0.intersection(factor_range*m0))
    out1 = out1.intersection(arb(0, (r*out0.upper()).upper()))
    return checked((out0, out1))


def enclose_two(node, h, lo, hi, cache):
    h, lo, hi = map(F, (h, lo, hi))
    if h <= 0 or hi <= lo:
        raise ValueError("positive mesh and nonempty cube required")
    if not node.children:
        key = (node.potential, h, lo, hi)
        if key not in cache:
            count = (hi-lo)/h
            if count.denominator != 1:
                raise ValueError("mesh and bounds must align")
            values = [leaf_two(node.potential, lo+k*h, lo+(k+1)*h,
                               lo+(F(k)+F(1,2))*h) for k in range(int(count))]
            cache[key] = TwoMessage([list(row) for row in zip(*values)], 1, lo)
        return cache[key]
    groups = [enclose_two(child, h, lo, hi, cache) for child in node.children]
    while len(groups) > 1:
        groups = [merge_two(groups[i], groups[i+1]) if i+1 < len(groups) else groups[i]
                  for i in range(0, len(groups), 2)]
    group = groups[0]
    radius = group.leaves*h/2
    result = [[], []]
    for index, moments in enumerate(zip(*group.moments)):
        center = group.origin+index*h+radius
        values = factor_two(moments, node.potential, center, radius)
        for row, value in zip(result, values):
            row.append(value)
    return TwoMessage(result, group.leaves, group.origin)


def partition_two(node, h, hi):
    h, hi = F(h), F(hi)
    message = enclose_two(node, h, F(0), hi, {})
    mass = sum(message.moments[0], arb(0))
    lower = mass.lower().max(0)
    tail = tail_envelope(node, hi)
    upper = (mass.upper()+tail).upper()
    width = (upper.log()-lower.log()).upper() if lower > 0 else arb.pos_inf()
    return lower, upper, tail, width, len(message.moments[0])
