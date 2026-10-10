"""Equivalent factor enclosures using exact quadratic recurrences on a mesh.

This only accelerates evaluation of the existing cell bounds. It does not
change their first-order discretization error, tree model, cube, or tail.
Arb encloses every exponential and recurrence multiplication.
"""
from fractions import Fraction as F

from flint import arb

from integration_bounds import (Message, ball, convolve_bounds, enclose,
                                nonnegative_lower, tail_envelope)


def ceil_fraction(x):
    return -((-x.numerator)//x.denominator)


def factor_bounds_on_grid(potential, h, origin, count, span):
    """For k=0..count-1, bound exp(-V) on [origin+kh,origin+(k+span)h]."""
    if h <= 0 or count <= 0 or span < 0:
        raise ValueError('positive mesh/count and nonnegative span required')
    stop=count+span
    # A branch starts at the first grid point on or above each target.
    cuts=sorted({0,stop,*[ceil_fraction((t-origin)/h) for t in potential.knots
                         if 0 < ceil_fraction((t-origin)/h) < stop]})
    factors=[]
    for start,end in zip(cuts,cuts[1:]):
        x=origin+start*h
        a,b,c=potential.coefficients(x)
        first=(-(ball(a*x*x-2*b*x+c))).exp()
        ratio=(-(ball(a*(2*x*h+h*h)-2*b*h))).exp()
        ratio_step=(-(ball(2*a*h*h))).exp()
        # Exact identity for exp of a quadratic; ball arithmetic propagates error.
        for k in range(start,end):
            factors.append(first)
            if k+1<end:
                first=first*ratio
                ratio=ratio*ratio_step
    assert len(factors)==stop
    q=(potential.minimizer-origin)/h
    # Intervals with k <= q <= k+span contain a global maximum of the factor.
    contains_start=ceil_fraction(q-span)
    contains_end=q.numerator//q.denominator
    maximum=(-potential.value(potential.minimizer)).exp().upper()
    lowers=[]
    uppers=[]
    for k in range(count):
        left,right=factors[k],factors[k+span]
        if k < contains_start:  # Interval entirely below the minimizer.
            lower,upper=nonnegative_lower(left),right.upper()
        elif k > contains_end:  # Interval entirely above the minimizer.
            lower,upper=nonnegative_lower(right),left.upper()
        else:
            lower,upper=nonnegative_lower(left.min(right)),maximum
        lowers.append(lower)
        uppers.append(upper)
    return lowers,uppers


def enclose_grid(node,h,lo,hi,cache):
    if not node.children:
        return enclose(node,h,lo,hi,cache)
    children=[enclose_grid(c,h,lo,hi,cache) for c in node.children]
    while len(children)>1:
        children=[convolve_bounds(children[i],children[i+1]) if i+1<len(children)
                  else children[i] for i in range(0,len(children),2)]
    combined=children[0]
    fl,fu=factor_bounds_on_grid(node.potential,h,combined.origin,
                               len(combined.lower),combined.n)
    return Message([nonnegative_lower(a*b) for a,b in zip(combined.lower,fl)],
                   [(a*b).upper() for a,b in zip(combined.upper,fu)],
                   combined.n,combined.origin)


def partition_enclosure_grid(node,h,hi,cache=None):
    msg=enclose_grid(node,h,F(0),hi,{} if cache is None else cache)
    tail=tail_envelope(node,hi)
    low=sum(msg.lower,arb(0)).lower().max(0)
    up=(sum(msg.upper,arb(0))+tail).upper()
    if not low>0:
        raise ArithmeticError('positive normalization lower bound unavailable')
    width=(up.log()-low.log()).upper()
    return low,up,tail,width,len(msg.lower)
