import Mathlib.Tactic

/- Structural containment for the finite moment operations used at a fixed mesh.
The indexed trace records semantic real moments. Leaf integration and the
factor's analytic remainder are explicit contracts, expressed through
Gaussian integral assumptions and MomentEnclosure. No Arb implementation is formalized.
The concrete arithmetic below uses exact real interval endpoints. -/
namespace MomentCircuit
noncomputable section

structure Bounds where
  lo : ℝ
  hi : ℝ

def Bounds.Contains (b : Bounds) (x : ℝ) : Prop := b.lo ≤ x ∧ x ≤ b.hi
def Bounds.add (a b : Bounds) : Bounds := ⟨a.lo + b.lo, a.hi + b.hi⟩
def Bounds.scale (c : ℝ) (a : Bounds) : Bounds :=
  ⟨min (c * a.lo) (c * a.hi), max (c * a.lo) (c * a.hi)⟩
def Bounds.mul (a b : Bounds) : Bounds :=
  ⟨min (min (a.lo*b.lo) (a.hi*b.lo)) (min (a.lo*b.hi) (a.hi*b.hi)),
   max (max (a.lo*b.lo) (a.hi*b.lo)) (max (a.lo*b.hi) (a.hi*b.hi))⟩
def Bounds.inflate (a : Bounds) (r : ℝ) : Bounds := ⟨a.lo-r, a.hi+r⟩
def Bounds.meet (a b : Bounds) : Bounds := ⟨max a.lo b.lo, min a.hi b.hi⟩

theorem add_contains {a b : Bounds} {x y : ℝ}
    (hx : a.Contains x) (hy : b.Contains y) : (a.add b).Contains (x+y) :=
  ⟨add_le_add hx.1 hy.1, add_le_add hx.2 hy.2⟩

theorem scale_contains {a : Bounds} {x c : ℝ} (hx : a.Contains x) :
    (a.scale c).Contains (c*x) := by
  unfold Bounds.Contains Bounds.scale
  by_cases hc : 0 ≤ c
  · exact ⟨(min_le_left _ _).trans (mul_le_mul_of_nonneg_left hx.1 hc),
      (mul_le_mul_of_nonneg_left hx.2 hc).trans (le_max_right _ _)⟩
  · have hc' : c ≤ 0 := le_of_lt (lt_of_not_ge hc)
    exact ⟨(min_le_right _ _).trans (mul_le_mul_of_nonpos_left hx.2 hc'),
      (mul_le_mul_of_nonpos_left hx.1 hc').trans (le_max_left _ _)⟩

theorem mul_contains {a b : Bounds} {x y : ℝ}
    (hx : a.Contains x) (hy : b.Contains y) : (a.mul b).Contains (x*y) := by
  have h0 := scale_contains (c := b.lo) hx
  have h1 := scale_contains (c := b.hi) hx
  have h := scale_contains (c := x) hy
  simp only [Bounds.Contains, Bounds.scale, mul_comm b.lo, mul_comm b.hi] at h0 h1 h
  unfold Bounds.Contains Bounds.mul
  constructor
  · exact (min_le_min h0.1 h1.1).trans h.1
  · exact h.2.trans (max_le_max h0.2 h1.2)

theorem inflate_contains {a : Bounds} {x y r : ℝ}
    (hx : a.Contains x) (herr : |y-x| ≤ r) : (a.inflate r).Contains y := by
  have hh := abs_le.mp herr
  constructor <;> dsimp [Bounds.inflate] <;> linarith [hx.1, hx.2, hh.1, hh.2]

theorem meet_contains {a b : Bounds} {x : ℝ}
    (ha : a.Contains x) (hb : b.Contains x) : (a.meet b).Contains x :=
  ⟨max_le ha.1 hb.1, le_min ha.2 hb.2⟩

structure Moments where
  mass : ℝ
  first : ℝ
  second : ℝ

structure MomentBounds where
  mass : Bounds
  first : Bounds
  second : Bounds

def MomentBounds.Contains (b : MomentBounds) (m : Moments) : Prop :=
  b.mass.Contains m.mass ∧ b.first.Contains m.first ∧ b.second.Contains m.second

def Moments.add (a b : Moments) : Moments :=
  ⟨a.mass+b.mass, a.first+b.first, a.second+b.second⟩
def MomentBounds.add (a b : MomentBounds) : MomentBounds :=
  ⟨a.mass.add b.mass, a.first.add b.first, a.second.add b.second⟩
def Moments.product (a b : Moments) : Moments :=
  ⟨a.mass*b.mass, a.first*b.mass+a.mass*b.first,
   a.second*b.mass+2*(a.first*b.first)+a.mass*b.second⟩
def MomentBounds.product (a b : MomentBounds) : MomentBounds :=
  ⟨a.mass.mul b.mass, (a.first.mul b.mass).add (a.mass.mul b.first),
   ((a.second.mul b.mass).add ((a.first.mul b.first).scale 2)).add
     (a.mass.mul b.second)⟩
def MomentBounds.meet (a b : MomentBounds) : MomentBounds :=
  ⟨a.mass.meet b.mass, a.first.meet b.first, a.second.meet b.second⟩

theorem add_moments_contains {a b : Moments} {A B : MomentBounds}
    (ha : A.Contains a) (hb : B.Contains b) : (A.add B).Contains (a.add b) :=
  ⟨add_contains ha.1 hb.1, add_contains ha.2.1 hb.2.1,
   add_contains ha.2.2 hb.2.2⟩

theorem product_moments_contains {a b : Moments} {A B : MomentBounds}
    (ha : A.Contains a) (hb : B.Contains b) : (A.product B).Contains (a.product b) :=
  ⟨mul_contains ha.1 hb.1,
   add_contains (mul_contains ha.2.1 hb.1) (mul_contains ha.1 hb.2.1),
   add_contains (add_contains (mul_contains ha.2.2 hb.1)
     (scale_contains (mul_contains ha.2.1 hb.2.1))) (mul_contains ha.1 hb.2.2)⟩

theorem meet_moments_contains {a : Moments} {A B : MomentBounds}
    (ha : A.Contains a) (hb : B.Contains a) : (A.meet B).Contains a :=
  ⟨meet_contains ha.1 hb.1, meet_contains ha.2.1 hb.2.1,
   meet_contains ha.2.2 hb.2.2⟩

/-- The three analytic bounds obtained by integration of the local remainders.
The output is the true factor-weighted moment triple, not a chosen approximation. -/
structure FactorContract (a b : Moments) (g d K1 K2 r : ℝ) : Prop where
  K1_nonneg : 0 ≤ K1
  K2_nonneg : 0 ≤ K2
  radius_nonneg : 0 ≤ r
  mass_error : |b.mass-(g*a.mass+d*a.first)| ≤ K2/2*a.second
  first_error : |b.first-(g*a.first+d*a.second)| ≤ K2*r/2*a.second
  second_error : |b.second-g*a.second| ≤ K1*r*a.second

def MomentBounds.factor (a : MomentBounds) (g d : Bounds) (K1 K2 r : ℝ) : MomentBounds :=
  ⟨((g.mul a.mass).add (d.mul a.first)).inflate (K2/2*a.second.hi),
   ((g.mul a.first).add (d.mul a.second)).inflate (K2*r/2*a.second.hi),
   (g.mul a.second).inflate (K1*r*a.second.hi)⟩

theorem factor_moments_contains {a b : Moments} {A : MomentBounds}
    {g d K1 K2 r : ℝ} {G D : Bounds} (ha : A.Contains a)
    (hg : G.Contains g) (hd : D.Contains d) (hc : FactorContract a b g d K1 K2 r) :
    (A.factor G D K1 K2 r).Contains b := by
  have e0 := hc.mass_error.trans
    (mul_le_mul_of_nonneg_left ha.2.2.2 (div_nonneg hc.K2_nonneg (by norm_num)))
  have e1 := hc.first_error.trans (mul_le_mul_of_nonneg_left ha.2.2.2
    (div_nonneg (mul_nonneg hc.K2_nonneg hc.radius_nonneg) (by norm_num)))
  have e2 := hc.second_error.trans (mul_le_mul_of_nonneg_left ha.2.2.2
    (mul_nonneg hc.K1_nonneg hc.radius_nonneg))
  exact ⟨inflate_contains (add_contains (mul_contains hg ha.1) (mul_contains hd ha.2.1)) e0,
    inflate_contains (add_contains (mul_contains hg ha.2.1) (mul_contains hd ha.2.2)) e1,
    inflate_contains (mul_contains hg ha.2.2) e2⟩

/-- Finite computation graph written as a tree; sharing affects cost, not soundness.
Every intermediate aggregate bin and final summation is represented by these
operations. A leaf is an analytic input enclosure; an intersection needs an
independent true-moment bound. The constructor assumptions are numerical and
analytic contracts, not assumptions that the final answer is enclosed. -/
inductive Trace : Moments → Type where
  | leaf (m : Moments) (b : MomentBounds) (h : b.Contains m) : Trace m
  | add {a b : Moments} (ta : Trace a) (tb : Trace b) : Trace (a.add b)
  | product {a b : Moments} (ta : Trace a) (tb : Trace b) : Trace (a.product b)
  | factor {a : Moments} (ta : Trace a) (b : Moments)
      (g d K1 K2 r : ℝ) (G D : Bounds) (hg : G.Contains g) (hd : D.Contains d)
      (hc : FactorContract a b g d K1 K2 r) : Trace b
  | meet {a : Moments} (ta : Trace a) (b : MomentBounds) (h : b.Contains a) : Trace a

def Trace.enclosure {a : Moments} : Trace a → MomentBounds
  | .leaf _ b _ => b
  | .add ta tb => ta.enclosure.add tb.enclosure
  | .product ta tb => ta.enclosure.product tb.enclosure
  | .factor ta _ _ _ K1 K2 r G D _ _ _ => ta.enclosure.factor G D K1 K2 r
  | .meet ta b _ => ta.enclosure.meet b

theorem trace_contains {a : Moments} (t : Trace a) : t.enclosure.Contains a := by
  induction t with
  | leaf m b h => exact h
  | add ta tb ha hb => exact add_moments_contains ha hb
  | product ta tb ha hb => exact product_moments_contains ha hb
  | factor ta b g d K1 K2 r G D hg hd hc ha =>
      exact factor_moments_contains ha hg hd hc
  | meet ta b hb ha => exact meet_moments_contains ha hb

/-- Dividing two positive/nonnegative integral enclosures. No independence of
the numerator and denominator is assumed or needed. -/
theorem ratio_enclosure {n d nl nu dl du : ℝ}
    (hn0 : 0 ≤ n) (hn : nl ≤ n ∧ n ≤ nu) (hd : dl ≤ d ∧ d ≤ du)
    (hdl : 0 < dl) : nl/du ≤ n/d ∧ n/d ≤ nu/dl := by
  have hd0 : 0 < d := hdl.trans_le hd.1
  have hdu : 0 < du := hd0.trans_le hd.2
  have hnu : 0 ≤ nu := hn0.trans hn.2
  constructor
  · exact (div_le_div_of_nonneg_right hn.1 hdu.le).trans
      (div_le_div_of_nonneg_left hn0 hd0 hd.2)
  · exact (div_le_div_of_nonneg_right hn.2 hd0.le).trans
      (div_le_div_of_nonneg_left hnu hdl hd.1)

end
end MomentCircuit

#print axioms MomentCircuit.mul_contains
#print axioms MomentCircuit.factor_moments_contains
#print axioms MomentCircuit.trace_contains
#print axioms MomentCircuit.ratio_enclosure
