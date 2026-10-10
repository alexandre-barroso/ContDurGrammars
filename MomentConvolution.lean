import Mathlib.MeasureTheory.Integral.Prod
import Mathlib.MeasureTheory.Integral.Bochner.Set
import Mathlib.Tactic

/- Continuous product-measure components of the moment argument. These prove
standard identities and a tail-envelope implication, not the program, Gaussian
primitives, coordinate partition, or finite-precision arithmetic. -/
namespace MomentConvolution
open MeasureTheory

variable {α β : Type*} [MeasurableSpace α] [MeasurableSpace β]
variable {μ : Measure α} {ν : Measure β}

theorem first_moment_product [IsFiniteMeasure μ] [IsFiniteMeasure ν]
    {f : α → ℝ} {g : β → ℝ} (hf : Integrable f μ) (hg : Integrable g ν) :
    (∫ z, (f z.1 + g z.2) ∂μ.prod ν) =
      ν.real Set.univ * (∫ x, f x ∂μ) + μ.real Set.univ * (∫ y, g y ∂ν) := by
  rw [integral_add (hf.comp_fst ν) (hg.comp_snd μ), integral_fun_fst, integral_fun_snd]
  simp only [smul_eq_mul]

theorem second_moment_product [IsFiniteMeasure μ] [IsFiniteMeasure ν]
    {f : α → ℝ} {g : β → ℝ} (hf : Integrable f μ) (hg : Integrable g ν)
    (hf2 : Integrable (fun x => f x ^ 2) μ)
    (hg2 : Integrable (fun y => g y ^ 2) ν) :
    (∫ z, (f z.1 + g z.2) ^ 2 ∂μ.prod ν) =
      ν.real Set.univ * (∫ x, f x ^ 2 ∂μ) +
      2 * (∫ x, f x ∂μ) * (∫ y, g y ∂ν) +
      μ.real Set.univ * (∫ y, g y ^ 2 ∂ν) := by
  have hfun : (fun z : α × β => (f z.1 + g z.2) ^ 2) =
      (fun z => (f z.1 ^ 2 + 2 * (f z.1 * g z.2)) + g z.2 ^ 2) := by
    funext z
    ring
  have hsum : Integrable (fun z : α × β => f z.1 ^ 2 + 2 * (f z.1 * g z.2))
      (μ.prod ν) := (hf2.comp_fst ν).add ((hf.mul_prod hg).const_mul 2)
  rw [hfun, integral_add hsum
    (hg2.comp_snd μ), integral_add (hf2.comp_fst ν) ((hf.mul_prod hg).const_mul 2),
    integral_fun_fst (fun x => f x ^ 2), integral_const_mul, integral_prod_mul,
    integral_fun_snd (fun y => g y ^ 2)]
  simp only [smul_eq_mul]
  ring

theorem tail_dominated_by_envelope {f g : α → ℝ} {s : Set α}
    (hs : MeasurableSet s) (hf : Integrable f μ) (hg : Integrable g μ)
    (hfg : ∀ x, f x ≤ g x) :
    (∫ x, f x ∂μ) - (∫ x in s, f x ∂μ) ≤
      (∫ x, g x ∂μ) - (∫ x in s, g x ∂μ) := by
  have hcomp := integral_mono (hf.restrict (s := sᶜ)) (hg.restrict (s := sᶜ)) hfg
  have hfs := integral_add_compl hs hf
  have hgs := integral_add_compl hs hg
  linarith

theorem full_integral_enclosure {f g : α → ℝ} {s : Set α} {lower upper : ℝ}
    (hs : MeasurableSet s) (hf : Integrable f μ) (hg : Integrable g μ)
    (hpos : ∀ x, 0 ≤ f x) (hfg : ∀ x, f x ≤ g x)
    (hlower : lower ≤ ∫ x in s, f x ∂μ)
    (hupper : (∫ x in s, f x ∂μ) ≤ upper) :
    lower ≤ (∫ x, f x ∂μ) ∧
      (∫ x, f x ∂μ) ≤ upper + (∫ x, g x ∂μ) - (∫ x in s, g x ∂μ) := by
  have hcomp : 0 ≤ ∫ x in sᶜ, f x ∂μ := integral_nonneg hpos
  have hfs := integral_add_compl hs hf
  have htail := tail_dominated_by_envelope hs hf hg hfg
  constructor <;> linarith

end MomentConvolution

#print axioms MomentConvolution.first_moment_product
#print axioms MomentConvolution.second_moment_product
#print axioms MomentConvolution.tail_dominated_by_envelope
#print axioms MomentConvolution.full_integral_enclosure
