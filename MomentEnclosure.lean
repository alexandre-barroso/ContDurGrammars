import Mathlib.MeasureTheory.Integral.Bochner.Basic
import Mathlib.Tactic

/- Local continuous-integral proof obligations.
Pointwise
Taylor/Lipschitz bounds and integrability are explicit hypotheses. This file
does not yet prove the complete tree algorithm, its derivative-bound code,
truncation implementation, or its floating-point library. -/
namespace MomentEnclosure
open MeasureTheory
open scoped ENNReal

variable {α : Type*} [MeasurableSpace α] {μ : Measure α}

theorem integral_error_bound {f g e : α → ℝ}
    (hf : Integrable f μ) (hg : Integrable g μ) (he : Integrable e μ)
    (h : ∀ᵐ x ∂μ, |f x - g x| ≤ e x) :
    |(∫ x, f x ∂μ) - (∫ x, g x ∂μ)| ≤ ∫ x, e x ∂μ := by
  rw [← integral_sub hf hg]
  exact abs_integral_le_integral_abs.trans
    (integral_mono_ae (hf.sub hg).abs he h)

theorem mass_propagation [IsFiniteMeasure μ] {d f : α → ℝ} {a b K : ℝ}
    (hd : Integrable d μ) (hd2 : Integrable (fun x => d x ^ 2) μ)
    (hf : Integrable f μ)
    (hrem : ∀ᵐ x ∂μ, |f x - (a + b * d x)| ≤ K / 2 * d x ^ 2) :
    |(∫ x, f x ∂μ) - (a * μ.real Set.univ + b * (∫ x, d x ∂μ))| ≤
      K / 2 * (∫ x, d x ^ 2 ∂μ) := by
  have hlin : Integrable (fun x => a + b * d x) μ :=
    (integrable_const a).add (hd.const_mul b)
  have h := integral_error_bound hf hlin (hd2.const_mul (K / 2)) hrem
  simpa only [integral_add (integrable_const a) (hd.const_mul b),
    integral_const, smul_eq_mul, integral_const_mul, mul_comm (μ.real Set.univ) a] using h

theorem first_moment_propagation {d f : α → ℝ} {a b K r : ℝ}
    (hr : 0 ≤ r) (_hK : 0 ≤ K)
    (hd : Integrable d μ) (hd2 : Integrable (fun x => d x ^ 2) μ)
    (hf1 : Integrable (fun x => d x * f x) μ)
    (hsupp : ∀ᵐ x ∂μ, |d x| ≤ r)
    (hrem : ∀ᵐ x ∂μ, |f x - (a + b * d x)| ≤ K / 2 * d x ^ 2) :
    |(∫ x, d x * f x ∂μ) -
      (a * (∫ x, d x ∂μ) + b * (∫ x, d x ^ 2 ∂μ))| ≤
      K * r / 2 * (∫ x, d x ^ 2 ∂μ) := by
  have hlin : Integrable (fun x => a * d x + b * d x ^ 2) μ :=
    (hd.const_mul a).add (hd2.const_mul b)
  have hpoint : ∀ᵐ x ∂μ,
      |d x * f x - (a * d x + b * d x ^ 2)| ≤ K * r / 2 * d x ^ 2 := by
    filter_upwards [hsupp, hrem] with x hx ht
    calc
      |d x * f x - (a * d x + b * d x ^ 2)| =
          |d x| * |f x - (a + b * d x)| := by
        rw [← abs_mul]
        congr 1
        ring
      _ ≤ r * (K / 2 * d x ^ 2) :=
        mul_le_mul hx ht (abs_nonneg _) hr
      _ = K * r / 2 * d x ^ 2 := by ring
  have h := integral_error_bound hf1 hlin (hd2.const_mul (K * r / 2)) hpoint
  simpa only [integral_add (hd.const_mul a) (hd2.const_mul b),
    integral_const_mul] using h

theorem second_moment_propagation {d f : α → ℝ} {a K r : ℝ}
    (hK : 0 ≤ K) (hd2 : Integrable (fun x => d x ^ 2) μ)
    (hf2 : Integrable (fun x => d x ^ 2 * f x) μ)
    (hsupp : ∀ᵐ x ∂μ, |d x| ≤ r)
    (hlip : ∀ᵐ x ∂μ, |f x - a| ≤ K * |d x|) :
    |(∫ x, d x ^ 2 * f x ∂μ) - a * (∫ x, d x ^ 2 ∂μ)| ≤
      K * r * (∫ x, d x ^ 2 ∂μ) := by
  have hpoint : ∀ᵐ x ∂μ,
      |d x ^ 2 * f x - a * d x ^ 2| ≤ K * r * d x ^ 2 := by
    filter_upwards [hsupp, hlip] with x hx ht
    calc
      |d x ^ 2 * f x - a * d x ^ 2| = d x ^ 2 * |f x - a| := by
        calc
          _ = |d x ^ 2 * (f x - a)| := by congr 1; ring
          _ = _ := by rw [abs_mul, abs_of_nonneg (sq_nonneg (d x))]
      _ ≤ d x ^ 2 * (K * r) :=
        mul_le_mul_of_nonneg_left (ht.trans (mul_le_mul_of_nonneg_left hx hK))
          (sq_nonneg _)
      _ = K * r * d x ^ 2 := by ring
  have h := integral_error_bound hf2 (hd2.const_mul a)
    (hd2.const_mul (K * r)) hpoint
  simpa only [integral_const_mul] using h

/-- Support bounds justify the first-moment intersection without treating its
signed integrand as nonnegative. -/
theorem first_support_bound [IsFiniteMeasure μ] {d : α → ℝ} {r : ℝ}
    (hd : Integrable d μ) (hsupp : ∀ᵐ x ∂μ, |d x| ≤ r) :
    |∫ x, d x ∂μ| ≤ r * μ.real Set.univ := by
  have h := abs_integral_le_integral_abs.trans
    (integral_mono_ae hd.abs (integrable_const r) hsupp)
  simpa only [integral_const, smul_eq_mul, mul_comm (μ.real Set.univ) r] using h

theorem second_support_bounds [IsFiniteMeasure μ] {d : α → ℝ} {r : ℝ}
    (hr : 0 ≤ r) (hd2 : Integrable (fun x => d x ^ 2) μ)
    (hsupp : ∀ᵐ x ∂μ, |d x| ≤ r) :
    0 ≤ (∫ x, d x ^ 2 ∂μ) ∧ (∫ x, d x ^ 2 ∂μ) ≤ r^2 * μ.real Set.univ := by
  have hp : ∀ᵐ x ∂μ, d x ^ 2 ≤ r^2 := by
    filter_upwards [hsupp] with x hx
    simpa only [sq_abs] using (sq_le_sq₀ (abs_nonneg (d x)) hr).mpr hx
  constructor
  · exact integral_nonneg (fun _ => sq_nonneg _)
  · have h := integral_mono_ae hd2 (integrable_const (r^2)) hp
    simpa only [integral_const, smul_eq_mul, mul_comm (μ.real Set.univ) (r^2)] using h

theorem factor_range_mass [IsFiniteMeasure μ] {f : α → ℝ} {l u : ℝ}
    (hf : Integrable f μ) (hl : ∀ᵐ x ∂μ, l ≤ f x) (hu : ∀ᵐ x ∂μ, f x ≤ u) :
    l * μ.real Set.univ ≤ (∫ x, f x ∂μ) ∧
      (∫ x, f x ∂μ) ≤ u * μ.real Set.univ := by
  constructor
  · have h := integral_mono_ae (integrable_const l) hf hl
    simpa only [integral_const, smul_eq_mul, mul_comm (μ.real Set.univ) l] using h
  · have h := integral_mono_ae hf (integrable_const u) hu
    simpa only [integral_const, smul_eq_mul, mul_comm (μ.real Set.univ) u] using h

theorem factor_range_second {d f : α → ℝ} {l u : ℝ}
    (hd2 : Integrable (fun x => d x ^ 2) μ)
    (hf2 : Integrable (fun x => d x ^ 2 * f x) μ)
    (hl : ∀ᵐ x ∂μ, l ≤ f x) (hu : ∀ᵐ x ∂μ, f x ≤ u) :
    l * (∫ x, d x ^ 2 ∂μ) ≤ (∫ x, d x ^ 2 * f x ∂μ) ∧
      (∫ x, d x ^ 2 * f x ∂μ) ≤ u * (∫ x, d x ^ 2 ∂μ) := by
  have hl' : ∀ᵐ x ∂μ, l * d x ^ 2 ≤ d x ^ 2 * f x := by
    filter_upwards [hl] with x hx
    nlinarith [mul_le_mul_of_nonneg_left hx (sq_nonneg (d x))]
  have hu' : ∀ᵐ x ∂μ, d x ^ 2 * f x ≤ u * d x ^ 2 := by
    filter_upwards [hu] with x hx
    nlinarith [mul_le_mul_of_nonneg_left hx (sq_nonneg (d x))]
  constructor
  · simpa only [integral_const_mul] using integral_mono_ae (hd2.const_mul l) hf2 hl'
  · simpa only [integral_const_mul] using integral_mono_ae hf2 (hd2.const_mul u) hu'

end MomentEnclosure

#print axioms MomentEnclosure.integral_error_bound
#print axioms MomentEnclosure.mass_propagation
#print axioms MomentEnclosure.first_moment_propagation
#print axioms MomentEnclosure.second_moment_propagation
#print axioms MomentEnclosure.first_support_bound
#print axioms MomentEnclosure.second_support_bounds
#print axioms MomentEnclosure.factor_range_mass
#print axioms MomentEnclosure.factor_range_second
