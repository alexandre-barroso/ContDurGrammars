"""Range enclosures for continuous duration models.

Parameters are exact rationals. Arb evaluates knot-split Gaussian integrals,
nonlinear factor bounds and the unary exterior envelope. See METHOD.md.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction as F
from functools import cached_property
import itertools
import json
from pathlib import Path
import platform
import time

import flint
from flint import arb, arb_poly, ctx

ROOT = Path(__file__).resolve().parent


def ball(q):
    q = F(q)
    return arb(q.numerator) / q.denominator


def nonnegative_lower(x):
    return x.lower().max(0)


def fraction_string(x):
    """Lossless serialization of an exact Arb endpoint."""
    return str(x.fmpq())


@dataclass(frozen=True)
class Potential:
    # Terms are (target, above weight, below weight), in consistent units.
    terms: tuple[tuple[F, F, F], ...] = ()

    def __post_init__(self):
        # Canonicalize before ANY coefficient arithmetic. In particular,
        # Fraction(float) preserves the exact supplied binary value; arithmetic
        # on a float followed by conversion to Arb would lose the enclosure.
        terms = tuple(tuple(F(x) for x in term) for term in self.terms)
        if any(len(term) != 3 for term in terms):
            raise ValueError('each potential term needs target, above, below')
        object.__setattr__(self, 'terms', terms)
        if any(a < 0 or b < 0 for _, a, b in self.terms):
            raise ValueError("nonnegative penalties are required")

    @cached_property
    def knots(self):
        return tuple(sorted({t for t, _, _ in self.terms}))

    def coefficients(self, x):
        x = F(x)
        aa = bb = cc = F(0)
        for t, above, below in self.terms:
            w = above if x >= t else below
            aa += w
            bb += w*t
            cc += w*t*t
        return aa, bb, cc  # V(x) = A*x² - 2*B*x + C

    def exact_value(self, x):
        x = F(x)
        return sum(((a if x >= t else b)*(x-t)**2
                    for t, a, b in self.terms), F(0))

    def value(self, x):
        x = F(x)
        # The branch is chosen with exact rational inputs.
        return sum((ball(a if x >= t else b)*ball(x-t)**2
                    for t, a, b in self.terms), arb(0))

    @cached_property
    def minimizer(self):
        if not self.terms:
            return F(0)
        candidates = list(self.knots)
        edges = [None, *self.knots, None]
        for l, r in zip(edges, edges[1:]):
            sample = (r-1 if l is None else l+1 if r is None else (l+r)/2)
            a, b, _ = self.coefficients(sample)
            if a:
                x = b/a
                if (l is None or x >= l) and (r is None or x <= r):
                    candidates.append(x)
        return min(candidates, key=self.exact_value)

    def factor_bounds(self, lo, hi):
        """Convex V: minimum at its clamped minimizer, maximum at endpoints."""
        lo, hi = F(lo), F(hi)
        if lo > hi:
            raise ValueError("reversed interval")
        x = min(hi, max(lo, self.minimizer))
        maximum = self.value(lo).max(self.value(hi))
        minimum = self.value(x)
        return nonnegative_lower((-maximum).exp()), (-minimum).exp().upper()

    def gaussian_piece(self, lo, hi, coefficients):
        """Integral of exp(-Ax²+2Bx-C); hi=None means positive infinity."""
        lo, hi = F(lo), None if hi is None else F(hi)
        a, b, c = map(F, coefficients)
        if not a:
            if hi is None:
                raise ValueError("leaf envelope is not integrable")
            if b:
                raise ValueError("unexpected linear-only potential")
            return ball(hi-lo)*(-ball(c)).exp()
        aa = ball(a)
        sq = aa.sqrt()
        mean = b/a
        zlo = sq*ball(lo-mean)
        scale = (ball(b*b/a-c)).exp()*arb.pi().sqrt()/(2*sq)
        if hi is None:
            return scale*zlo.erfc()
        zhi = sq*ball(hi-mean)
        # Tail differences avoid cancellation of two values near +/-1.
        if lo >= mean:
            difference = zlo.erfc()-zhi.erfc()
        elif hi <= mean:
            difference = (-zhi).erfc()-(-zlo).erfc()
        else:
            difference = zhi.erf()-zlo.erf()
        return scale*difference

    def integral(self, lo, hi=None):
        lo, hi = F(lo), None if hi is None else F(hi)
        if hi is not None and hi < lo:
            raise ValueError("reversed limits")
        cuts = [lo, *(x for x in self.knots if x > lo and (hi is None or x < hi)), hi]
        total = arb(0)
        for left, right in zip(cuts, cuts[1:]):
            if right == left:
                continue
            sample = left+1 if right is None else (left+right)/2
            total += self.gaussian_piece(left, right, self.coefficients(sample))
        if not total.is_finite():
            raise ArithmeticError("nonfinite primitive")
        return total


@dataclass(frozen=True)
class Node:
    potential: Potential
    children: tuple[Node, ...] = ()
    name: str = ""

    @cached_property
    def leaves(self):
        if not self.children:
            return (self,)
        return tuple(itertools.chain.from_iterable(c.leaves for c in self.children))


@dataclass
class Message:
    lower: list
    upper: list
    n: int
    origin: F


def convolve_bounds(left, right):
    low = arb_poly(left.lower)*arb_poly(right.lower)
    up = arb_poly(left.upper)*arb_poly(right.upper)
    ncoeff = len(left.lower)+len(right.lower)-1
    # Even zero trailing coefficients count as cells in this geometry.
    lows = [nonnegative_lower(low[k]) for k in range(ncoeff)]
    ups = [up[k].upper() for k in range(ncoeff)]
    if any(x < 0 for x in ups):
        raise ArithmeticError("negative upper coefficient")
    return Message(lows, ups, left.n+right.n, left.origin+right.origin)


def enclose(node, h, lo, hi, cache):
    if not node.children:
        key = (node.potential, h, lo, hi)
        if key not in cache:
            count = (hi-lo)/h
            if count.denominator != 1:
                raise ValueError("bounds must align with mesh")
            masses = [node.potential.integral(lo+k*h, lo+(k+1)*h)
                      for k in range(int(count))]
            cache[key] = Message([nonnegative_lower(x) for x in masses],
                                 [x.upper() for x in masses], 1, lo)
        return cache[key]
    children = [enclose(c, h, lo, hi, cache) for c in node.children]
    # Balanced arithmetic reduction; it does not alter grammatical grouping.
    while len(children) > 1:
        children = [convolve_bounds(children[i], children[i+1])
                    if i+1 < len(children) else children[i]
                    for i in range(0, len(children), 2)]
    combined = children[0]
    lower, upper = [], []
    for k, (a, b) in enumerate(zip(combined.lower, combined.upper)):
        start = combined.origin+k*h
        lfactor, ufactor = node.potential.factor_bounds(start, start+combined.n*h)
        lower.append(nonnegative_lower(a*lfactor))
        upper.append((b*ufactor).upper())
    return Message(lower, upper, combined.n, combined.origin)


def tail_envelope(node, hi):
    """Outside [0,hi]^n: omit all nonnegative internal penalties."""
    full = inside = arb(1)
    for leaf in node.leaves:
        full *= leaf.potential.integral(F(0))
        inside *= leaf.potential.integral(F(0), hi)
    return (full-inside).upper().max(0)


def partition_enclosure(node, h, hi, cache=None):
    h, hi = F(h), F(hi)
    msg = enclose(node, h, F(0), hi, {} if cache is None else cache)
    tail = tail_envelope(node, hi)
    low = sum(msg.lower, arb(0)).lower().max(0)
    up = (sum(msg.upper, arb(0))+tail).upper()
    if not low > 0:
        raise ArithmeticError("positive normalization lower bound unavailable")
    width = (up.log()-low.log()).upper()
    return low, up, tail, width, len(msg.lower)


def hayes_models():
    # All model constants and their attribution are in the adjacent parameter file.
    p = json.loads(Path(__file__).with_name("published_timing_parameters.json").read_text())
    # Rational strings specify the input binary64 values exactly.
    w = [F(x)/p["harmony_divisor"] for x in p["weights_in_column_order_AU_to_BB"]]
    t = {k:F(v) for k,v in p["targets"].items()}
    structures = {"HH0H":((2,2),(2,2)), "HHLH":((2,2,1),(1,2,1)),
                  "HLLH":((2,1,1),(2,1,1)), "LHLH":((1,2,1),(1,2,1)),
                  "LLHH":((1,1,2),(1,1,2))}
    models = {}
    for name,(moras,grids) in structures.items():
        leaves = []
        for i,(m,g) in enumerate(zip(moras,grids)):
            pot = Potential(((m*t["mora"],w[0]/m,w[1]/m),
                             (t["syllable"],w[2],w[3]),
                             (g*t["grid"],w[4]/g,w[5]/g)))
            leaves.append(Node(pot,name=f"{name}_{i}"))
        root = Potential(((t["hemimetron"],w[6],w[7]),))
        models[name] = Node(root,tuple(leaves),name)
    return models


def benchmark(output, max_power):
    ctx.prec = 128
    start = datetime.now(timezone.utc).isoformat()
    result = {"started_at":start, "method":"leaf-exact integrals and sum-cell factor enclosures", "status":"NUMERICAL_FEASIBILITY_NOT_NOVELTY_OR_FORMAL_PROOF", "arithmetic":{"python_flint":flint.__version__,"precision_bits":ctx.prec,
                            "python":platform.python_version(),"machine":platform.machine()}, "models":{}, "acoustic_data_or_counts_used":False}
    output.parent.mkdir(parents=True,exist_ok=True)
    def save():
        output.write_text(json.dumps(result,indent=2)+"\n")
    save()
    for name,node in hayes_models().items():
        rows=[]
        cumulative_seconds=0.0
        result["models"][name]={"dimension":len(node.leaves),"attempts":rows}
        for j in range(max_power+1):
            h=F(20,2**j)
            before=time.perf_counter()
            low,up,tail,width,ncoeff=partition_enclosure(node,h,F(800))
            elapsed=time.perf_counter()-before
            cumulative_seconds += elapsed
            row={"mesh_ms":str(h),"seconds":elapsed,"coefficients_at_root":ncoeff,
                 "cumulative_seconds":cumulative_seconds,
                 "Z_lower_exact_rational":fraction_string(low),
                 "Z_upper_exact_rational":fraction_string(up),
                 "tail_upper_exact_rational":fraction_string(tail),
                 "log_width_upper_exact_rational":fraction_string(width),
                 "log_width_approx":float(width),"Z_mid_approx":float((low+up)/2),
                 "target_proved_by_arithmetic":bool(width <= ball(F(1,1000)))}
            rows.append(row)
            print(json.dumps({"model":name,**{k:row[k] for k in ("mesh_ms","seconds","log_width_approx","target_proved_by_arithmetic")}}),flush=True)
            save()
            if row["target_proved_by_arithmetic"] or cumulative_seconds > 60:
                break
        result["models"][name]["feasibility_target_met"]=bool(rows[-1]["target_proved_by_arithmetic"] and cumulative_seconds<=60)
        save()
    result["completed_at"]=datetime.now(timezone.utc).isoformat()
    save()


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--max-power",type=int,default=11)
    args=parser.parse_args()
    benchmark(args.output,args.max_power)
