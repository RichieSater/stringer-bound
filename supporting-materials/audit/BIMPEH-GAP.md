# Reassessment of Bimpeh's finite-sample certification of the Stringer bound

**Claim under audit.** Bimpeh (PhD thesis, Dublin City University, 2008,
ch. 5, [full text](https://doras.dcu.ie/600/1/YawThesis.PDF)) derives, for
the Stringer bound with binomial (Clopper–Pearson) factors, the coverage
lower bound

> CP ≥ P̄ₙ := P( U_{i:n} ≤ pₙ(i; 1−α), 1 ≤ i ≤ n )   (his eq. 5.16)

for taint distributions F on [0,1], computes P̄ₙ by Bolshev's recursion,
and concludes that the conservatism conjecture holds for n ≤ 11 at
α = 0.05 (n ≤ 10 at 0.1; …; n ≤ 7 at 0.5).

**Verdict: inequality (5.16) fails for both atomic and continuous F, and
the corrected containment event does not recover the claimed finite-sample
range.** P̄ₙ itself is correctly computed (this
repository reproduces his Table 5.1 to all printed digits;
`bolshev.py`), but it is not a lower bound on coverage.

## The error

Bimpeh's band is F̂ₙ,L(t) = qₙ(i−1; 1−α) on [t_{i−1:n}, t_{i:n}) (his
5.12), with qₙ(i) the lower CP limit for i successes, qₙ(0) = 0. His
chain (thesis p. 78) is

> CP ≥ P[F(t) ≥ F̂ₙ,L(t) ∀t] = P[∩ᵢ F(t_{i:n}) ≥ F̂ₙ,L(t_{i:n})]
> = P(U_{i:n} ≥ qₙ(i−1), 1 ≤ i ≤ n) = P̄ₙ.

The first equality is fine: since F is nondecreasing and the band is a
right-continuous step function, pointwise domination is equivalent to
domination at the jump points, where the band's value is
F̂ₙ,L(t_{i:n}) = **qₙ(i)** — the band on the interval *to the right* of
t_{i:n}. The next step substitutes **qₙ(i−1)** — the band's left limit —
silently weakening every constraint by one index and dropping the top
constraint entirely. The correct containment event is

    F(t_{i:n}) ≥ qₙ(i)  for i = 1, …, n,

whose i = n constraint is F(t_{n:n}) ≥ qₙ(n) = α^{1/n}.

## The corrected containment event

For continuous F (where F(t_{i:n}) = U_{i:n} is legitimate), the correct
containment probability satisfies

    P( U_{i:n} ≥ qₙ(i) ∀i ) ≤ P( U_{n:n} ≥ α^{1/n} ) = 1 − α,

since P(U_{n:n} < x) = xⁿ. So the pointwise-band route yields at most
1−α. Thus the corrected pointwise event does not supply the claimed
strict numerical cushion above nominal. If the conjecture is true more
generally, it may be because the *integral* inequality
∫F̂ₙ,L ≤ ∫F holds on a strictly larger event than pointwise
domination.

## Counterexamples to (5.16), fully by hand

Let α = 0.05 and let F put mass 0.99 uniformly on (0.99, 1) and mass
0.01 uniformly on (0, 0.01) — a continuous distribution. Then
μ = 0.99·0.995 + 0.01·0.005 = 0.9851.

- **n = 1**: p₁(0) = 0.95, p₁(1) = 1, so SB = 0.95 + 0.05·T₁ and P̄₁ =
  P(U ≤ p₁(1)) = 1. SB < μ ⟺ T₁ < (0.9851−0.95)/0.05 = 0.702, which
  happens exactly when T₁ falls in the low cluster:
  **CP = 0.99 < 1 = P̄₁.**
- **n = 2**: p₂(0) = 1−√0.05 ≈ 0.776393, p₂(1) = √0.95 ≈ 0.974679,
  p₂(2) = 1; P̄₂ = 2√0.95 − 0.95 ≈ 0.999359. With S ≥ T the two sampled
  taints, SB = p₂(0) + (p₂(1)−p₂(0))S + (1−p₂(1))T. If both taints are
  high (probability 0.99²), SB ≥ 0.776 + 0.198·0.99 + 0.025·0.99 >
  0.9851: covered. If S is high and T low (probability 2·0.99·0.01 =
  0.0198), SB ≤ 0.776 + 0.198·1 + 0.025·0.01 ≈ 0.9749 < 0.9851: fails.
  Both low: fails. **CP = 1 − 0.0199 = 0.9801 < 0.99936 = P̄₂.**

Both are confirmed by Monte Carlo (4·10⁶ samples, agreement to 4
decimals; `bimpeh_continuous_check.py`). An exact atomic counterexample
at n = 5 is also on record: P(T=1) = 1/2, P(T=1/10) = 1/4, P(T=0) = 1/4
has CP = 31/32 = 0.96875 < P̄₅ = 0.987458, in exact rational arithmetic
(`bolshev.py --crosscheck`), with the failure event being exactly the
1/32-probability event of drawing no taint-1 unit.

Note the conjecture itself survives in every counterexample (0.99,
0.9801, 0.96875 ≥ 0.95): what fails is Bimpeh's *certificate*, not the
Stringer bound.

## Consequences

1. Bimpeh's Table 5.1 frontier ("reliable for n ≤ 11 at α ≤ 0.05") does
   not certify the stated class of all F. His n ≤ 2 "analytic proof"
   establishes P̄₂ ≥ 1−α, a true statement about a quantity that is not
   the claimed coverage bound.
2. Apart from the \(n=2\) theorem proved in this repository, the established
   finite-sample cases include n = 1 (Bickel 1992),
   P(SB ≥ μ) ≥ (1−α)^{n+1} for n ≥ 2 "under certain conditions on F"
   (Bickel 1992, hypotheses unverified), the {0,1}-support case
   (de Jager–Pap–van Zuijlen 1997), and the single-nonzero-taint case
   (`two_point_lemma.py`, this repository).
3. Numerical searches in this repository over distributions with two and
   three nonzero taint values (allowing mass at zero) found no violation at
   α = 0.05 for the reported sample sizes. The smallest coverage found
   agreed with 1−α to search precision.
