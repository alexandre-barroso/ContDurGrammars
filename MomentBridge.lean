import MomentEnclosure
import MomentCircuit

/- This bridge discharges the circuit's analytic factor contract from actual
integrals. The partition/Fubini representation is justified separately in the
readable global proof; Gaussian and Arb software remain trusted contracts. -/
namespace MomentBridge
open MeasureTheory MomentCircuit
noncomputable section

variable {α : Type*} [MeasurableSpace α] {μ : Measure α} [IsFiniteMeasure μ]

def before (d : α → ℝ) (μ : Measure α) : Moments :=
  ⟨μ.real Set.univ, ∫ x, d x ∂μ, ∫ x, d x^2 ∂μ⟩
def after (d f : α → ℝ) (μ : Measure α) : Moments :=
  ⟨∫ x, f x ∂μ, ∫ x, d x*f x ∂μ, ∫ x, d x^2*f x ∂μ⟩

theorem factor_contract_from_integrals {d f : α → ℝ} {g b K1 K2 r : ℝ}
    (hK1 : 0 ≤ K1) (hK2 : 0 ≤ K2) (hr : 0 ≤ r)
    (hd : Integrable d μ) (hd2 : Integrable (fun x => d x^2) μ)
    (hf : Integrable f μ) (hf1 : Integrable (fun x => d x*f x) μ)
    (hf2 : Integrable (fun x => d x^2*f x) μ)
    (hsupp : ∀ᵐ x ∂μ, |d x| ≤ r)
    (hlip : ∀ᵐ x ∂μ, |f x-g| ≤ K1*|d x|)
    (hrem : ∀ᵐ x ∂μ, |f x-(g+b*d x)| ≤ K2/2*d x^2) :
    FactorContract (before d μ) (after d f μ) g b K1 K2 r := by
  exact ⟨hK1, hK2, hr,
    MomentEnclosure.mass_propagation hd hd2 hf hrem,
    MomentEnclosure.first_moment_propagation hr hK2 hd hd2 hf1 hsupp hrem,
    MomentEnclosure.second_moment_propagation hK1 hd2 hf2 hsupp hlip⟩

theorem enclosed_factor_from_integrals {d f : α → ℝ} {g b K1 K2 r : ℝ}
    (hK1 : 0 ≤ K1) (hK2 : 0 ≤ K2) (hr : 0 ≤ r)
    (hd : Integrable d μ) (hd2 : Integrable (fun x => d x^2) μ)
    (hf : Integrable f μ) (hf1 : Integrable (fun x => d x*f x) μ)
    (hf2 : Integrable (fun x => d x^2*f x) μ)
    (hsupp : ∀ᵐ x ∂μ, |d x| ≤ r)
    (hlip : ∀ᵐ x ∂μ, |f x-g| ≤ K1*|d x|)
    (hrem : ∀ᵐ x ∂μ, |f x-(g+b*d x)| ≤ K2/2*d x^2)
    (A : MomentBounds) (G B : Bounds) (ha : A.Contains (before d μ))
    (hg : G.Contains g) (hb : B.Contains b) :
    (A.factor G B K1 K2 r).Contains (after d f μ) := by
  exact factor_moments_contains ha hg hb
    (factor_contract_from_integrals hK1 hK2 hr hd hd2 hf hf1 hf2 hsupp hlip hrem)

end
end MomentBridge

#print axioms MomentBridge.factor_contract_from_integrals
#print axioms MomentBridge.enclosed_factor_from_integrals
