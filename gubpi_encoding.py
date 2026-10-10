"""SPCF encodings and output parsing for the probabilistic-program baseline.

Arithmetic and Gaussian-PDF encodings expose their finite comparison checks.
See GUBPI.md for execution and model-equivalence limits.
"""
from fractions import Fraction as F
import ast
import math
import re
import mpmath as mp
from integration_bounds import ball, fraction_string


def integer_expression(n):
    if n < 0:
        return f'(-{integer_expression(-n)})'
    if n <= 2**53:
        return str(n)
    q, r = divmod(n, 10**9)
    return f'({integer_expression(q)}*1000000000+{r})'

def rational_expression(x):
    x = F(x)
    if x.denominator == 1:
        return integer_expression(x.numerator)
    return f'({integer_expression(x.numerator)}/{integer_expression(x.denominator)})'

def proposal_parameters(node):
    result=[]
    for leaf in node.leaves:
        center=F(round(leaf.potential.minimizer/800*1024),1024)
        curvature=sum((min(a,b)*800**2 for _,a,b in leaf.potential.terms),F(0))
        if curvature<=0:
            raise ValueError('Gaussian proposal rule requires positive minimum unary curvature')
        sigma=F(1,128)
        while 1/(2*sigma**2)>curvature:
            sigma*=2
        result.append((center,sigma))
    return result

def arithmetic_program(node, mode, setting, proposal='uniform'):
    """Return SPCF and a transparent list of arithmetic expressions to audit."""
    lines = [f'# method {mode}', '# depth 2000', '# discretization -0.5 0.5 1',
             '# outputSplitProgress 0', '# outputCurrentPath 0', '# outputCurrentArea 0']
    if mode == 'boxes':
        lines += [f'# splits {setting}']
    else:
        lines += [f'# epsilonScore {setting}', f'# epsilonVar {setting}']
    proposals=proposal_parameters(node) if proposal=='gaussian' else None
    for i in range(len(node.leaves)):
        if proposals:
            mu,sigma=proposals[i]
            lines.append(f'let x{i} = sample normal({rational_expression(mu)},{rational_expression(sigma)}) in')
            lines += [f'if (x{i}<0) then (score(0);0) else',
                      f'if (x{i}>1) then (score(0);0) else']
        else:
            lines.append(f'let x{i} = sample uniform(0,1) in')
    leaves = {id(leaf): i for i, leaf in enumerate(node.leaves)}
    records = []

    def visit(v):
        for child in v.children:
            visit(child)
        indices = [leaves[id(leaf)] for leaf in v.leaves]
        s = '(' + '+'.join(f'x{i}' for i in indices) + ')'
        defs, costs = [], []
        group = len(records)
        for j, (target, above, below) in enumerate(v.potential.terms):
            d, a, pos, neg = (f'{prefix}{group}_{j}' for prefix in ('d', 'a', 'p', 'n'))
            defs += [(d, f'({s}-{rational_expression(target/800)})'),
                     (a, f'sqrt({d}*{d})'), (pos, f'(({d}+{a})/2)'),
                     (neg, f'((-{d}+{a})/2)')]
            costs.append(f'({rational_expression(above*800**2)}*{pos}*{pos}'
                         f'+{rational_expression(below*800**2)}*{neg}*{neg})')
        for name, expression in defs:
            lines.append(f'let {name} = {expression} in')
        cost = '(' + '+'.join(costs) + ')' if costs else '0'
        correction=None
        if proposals and not v.children:
            i=indices[0];mu,sigma=proposals[i]
            d=f'(x{i}-{rational_expression(mu)})'
            cost=f'({cost}-({d}*{d}/{rational_expression(2*sigma**2)}))'
            correction={'mean':str(mu),'standard_deviation':str(sigma)}
        lines.append(f'score(exp(-{cost}));')
        records.append({'indices': indices, 'definitions': defs, 'cost': cost,
                        'terms': [[str(x) for x in term] for term in v.potential.terms],
                        'proposal_correction':correction})
    visit(node)
    lines.append('0')
    return '\n'.join(lines) + '\n', records

def exact_expression(text, variables):
    """Independent Python-AST evaluator of emitted real arithmetic at rationals."""
    def ev(n):
        if isinstance(n, ast.Constant) and isinstance(n.value, int):
            assert abs(n.value) <= 2**53, 'inexact numeric literal'
            return F(n.value)
        if isinstance(n, ast.Name):
            return variables[n.id]
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
            return -ev(n.operand)
        if isinstance(n, ast.BinOp):
            a, b = ev(n.left), ev(n.right)
            if isinstance(n.op, ast.Add): return a+b
            if isinstance(n.op, ast.Sub): return a-b
            if isinstance(n.op, ast.Mult): return a*b
            if isinstance(n.op, ast.Div): return a/b
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == 'sqrt':
            q = ev(n.args[0]); assert q >= 0
            a, b = math.isqrt(q.numerator), math.isqrt(q.denominator)
            assert a*a == q.numerator and b*b == q.denominator
            return F(a,b)
        raise ValueError(ast.dump(n))
    return ev(ast.parse(text, mode='eval').body)

def check_arithmetic_encoding(node, records):
    count = 0
    # All factor knots and nearby points, rather than only a favorable mode.
    for record in records:
        terms = [[F(x) for x in t] for t in record['terms']]
        sums = {F(0), F(1), F(800*len(record['indices']))}
        for target, _, _ in terms:
            sums.update((target-1, target, target+1))
        for s in sums:
            variables = {f'x{i}': F(0) for i in range(len(node.leaves))}
            variables[f"x{record['indices'][0]}"] = s/800
            for name, expression in record['definitions']:
                variables[name] = exact_expression(expression, variables)
            actual = exact_expression(record['cost'], variables)
            expected = sum(((above if s >= target else below)*(s-target)**2
                            for target, above, below in terms), F(0))
            if record.get('proposal_correction'):
                q=record['proposal_correction'];mu=F(q['mean']);sigma=F(q['standard_deviation'])
                expected-=(s/800-mu)**2/(2*sigma*sigma)
            assert actual == expected, (actual, expected, s)
            count += 1
    return count

def parse_output(text, dimension, tail, proposal_constant=None):
    # GuBPI prints shortest round-trip decimals, not exact rational endpoints.
    m = re.search(r'^Normalizing constant: (\S+) <= Z <= (\S+)\s*$', text, re.M)
    if not m:
        return {'status':'NO_FINITE_NORMALIZER_OUTPUT'}
    a, b = float(m[1]), float(m[2])
    if not (math.isfinite(a) and math.isfinite(b) and 0 <= a <= b and b > 0):
        return {'status':'INVALID_OR_NONFINITE_NORMALIZER', 'printed_endpoints':list(m.groups())}
    # One binary64 step covers parsing/printing; does not certify the backend.
    lower = F(math.nextafter(a, -math.inf)) if a else F(0)
    upper = F(math.nextafter(b, math.inf))
    scale=ball(F(800)**dimension)
    if proposal_constant is not None:
        scale*=proposal_constant
    cube_lower=(ball(lower)*scale).lower().max(0)
    cube_upper=(ball(upper)*scale).upper()
    total_upper = (cube_upper+tail).upper()
    width = (total_upper.log()-cube_lower.log()).upper() if cube_lower>0 else None
    return {'status':'REPORTED_BOUNDS_PARSED', 'printed_endpoints':list(m.groups()),
            'cube_lower_from_printed_rational':fraction_string(cube_lower),
            'cube_upper_from_printed_rational':fraction_string(cube_upper),
            'tail_upper_exact_rational':fraction_string(tail),
            'full_domain_reported_log_width':float(width) if width is not None else None,
            'reported_target_met': bool(width is not None and width <= ball(F(1,1000))),
            'certificate_status':'PUBLISHED_FLOATING_IMPLEMENTATION_NOT_INDEPENDENTLY_CERTIFIED'}

def m(x):
    x=F(x);return mp.mpf(x.numerator)/x.denominator

def gaussian_branches(node):
    size=len(node.leaves)
    cuts=[F(0),*(t/800 for t in node.potential.knots if 0<t<size*800),F(size)]
    pieces=[]
    for l,r in zip(cuts,cuts[1:]):
        a,b,c=node.potential.coefficients((l+r)*400)
        a,b=a*800**2,b*800
        if a:
            mean=b/a;sd=1/mp.sqrt(2*m(a));amplitude=mp.exp(m(b*b/a-c))*mp.sqrt(mp.pi/m(a))
            pieces.append({'left':str(l),'right':str(r),'coefficients':[str(a),str(b),str(c)],
                           'mean':float(m(mean)),'sd':float(sd),'amplitude':float(amplitude)})
        else:
            if b:raise ValueError('linear-only branch unsupported')
            pieces.append({'left':str(l),'right':str(r),'coefficients':[str(a),str(b),str(c)],
                           'constant':float(mp.exp(-m(c)))})
    return pieces

def gaussian_program(node,epsilon):
    mp.mp.dps=90
    lines=['# method polytopes','# depth 2000','# discretization -0.5 0.5 1',
           '# outputSplitProgress 0','# outputCurrentPath 0','# outputCurrentArea 0',
           f'# epsilonScore {epsilon}',f'# epsilonVar {epsilon}']
    for i in range(len(node.leaves)):lines.append(f'let x{i} = sample uniform(0,1) in')
    indices={id(leaf):i for i,leaf in enumerate(node.leaves)};records=[]
    def visit(v):
        for child in v.children:visit(child)
        ids=[indices[id(leaf)] for leaf in v.leaves];s='('+ '+'.join(f'x{i}' for i in ids)+')'
        branches=gaussian_branches(v)
        def score(piece):
            if 'constant' in piece:return f'score({piece["constant"]:.17g})'
            return f'score(({piece["amplitude"]:.17g})*pdfnormal({piece["mean"]:.17g},{piece["sd"]:.17g},{s}))'
        expression=score(branches[-1])
        for piece in reversed(branches[:-1]):
            expression=f'if ({s} < {float(F(piece["right"])):.17g}) then ({score(piece)}) else ({expression})'
        lines.append('('+expression+');')
        records.append({'name':v.name,'indices':ids,'pieces':branches})
    visit(node);lines.append('0')
    return '\n'.join(lines)+'\n',records

def check_gaussian_encoding(records):
    count=0;max_relative=mp.mpf(0);max_guard=mp.mpf(0)
    for record in records:
        for p in record['pieces']:
            lo,hi=F(p['left']),F(p['right']);a,b,c=map(F,p['coefficients'])
            max_guard=max(max_guard,abs(m(F(float(hi)))-m(hi)))
            for x in [lo,(2*lo+hi)/3,(lo+hi)/2,(lo+2*hi)/3,hi]:
                exact=mp.exp(-m(a)*m(x)**2+2*m(b)*m(x)-m(c))
                if 'constant' in p:encoded=m(F(p['constant']))
                else:
                    sd,mu,amp=map(lambda k:m(F(p[k])),['sd','mean','amplitude'])
                    encoded=amp/(sd*mp.sqrt(2*mp.pi))*mp.exp(-(m(x)-mu)**2/(2*sd*sd))
                max_relative=max(max_relative,abs(encoded/exact-1));count+=1
    # A bounded encoding sanity check, not a uniform error proof.
    if max_relative>=mp.mpf('1e-10'):raise ArithmeticError('Gaussian encoding discrepancy exceeds sanity threshold')
    return {'points':count,'max_sampled_relative_factor_discrepancy':mp.nstr(max_relative,25),
            'max_guard_absolute_rounding':mp.nstr(max_guard,25),
            'scope':'finite high-precision checks of emitted float parameters; not a uniform model-equivalence certificate'}
