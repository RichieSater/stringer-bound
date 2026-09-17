# Novelty and priority verification

The priority assessment uses targeted literature searches, OpenAlex citation
tracing, Crossref metadata, and primary-source reading, including the full
Bimpeh thesis, the 2021 *Statistical Papers* MUS survey, the Dutch
statistical-auditing essay collection, and the 2024
Balakrishnan--Rychlik monograph. Initial date: 2026-08-07; searches were
extended on 2026-09-02 and 2026-09-08. The 2026-09-17 update checks the
external Gaffke theorem and related simplex-cap literature against primary
sources and records the unresolved active-gap implication for the proposed
five- and six-coordinate global theorems. The principal manuscript now
focuses on the four-coordinate theorem, its sharp n=3 Stringer consequence,
and the all-dimensional terminal-edge obstruction. The other results below
belong to the supporting material or the systematic-PPS companion.
Priority statements are qualified by the scope and date of the searches.

## Claim 1 — qualified priority statement for the n = 2 theorem

The literature through 2021 describes the finite-sample problem as open.
Bimpeh's thesis §5.3 claims the (n\le2) result through
P̄₂ = 2√(1−α) − (1−α) ≥ 1−α, but the supporting correction gives
continuous counterexamples to the containment inequality at (n=1,2).
Consequently, the cited argument does not establish the theorem. The searched
sources support qualified priority wording; the theorem
itself rests on the written wedge proof rather than on this literature
assessment.

## Claim 2 — the wedge inequality is new: HOLDS

- **Rychlik school** (sharp order-statistic/L-statistic bounds;
  Balakrishnan–Rychlik 2024 monograph): covers expectations and
  distribution functions of single order statistics — no sharp bound on
  P(L-statistic ≥ t) under a mean constraint exists there.
- **Closest relative**: Hoeffding & Shrikhande (1955) solved the
  equal-weights case (sharp bound for P(X₁+X₂ ≥ t), nonnegative iid,
  given mean); Łuczak–Mieczkowska–Šileikis (arXiv:1602.03547) extended
  to general k. The wedge inequality has unequal min/max weights except
  at the interior value α = 16/25. Both now cited.
- **Valid-mean-bound genre** (Anderson 1969; Gaffke 2005;
  Learned-Miller–Thomas; Phan et al. 2021; Bentkus–van Zuijlen 2003;
  Vlassis–Thomas 2026 proof of Gaffke's conjecture): none proves or
  implies the wedge inequality; Phan et al. cite the Stringer
  literature and list its coverage as unknown.
- **Trapdoor check (performed, negative)**: could n = 2 conservatism
  follow from a known-valid bound pointwise below SB?
  - Anderson's bound: no — SB fails to dominate it (gap −0.0253 at
    sample (0,1), α = 0.05).
  - Gaffke's bound: at α = 0.05, SB *does* pointwise dominate it
    (minimum gap 0, ties exactly on the two extremal families:
    equal-taint samples and samples with one taint = 1) — so at
    standard confidence an alternative route would exist *given* the
    2026 validity proof plus a proof of the (itself new) domination;
    but at α = 0.7 the domination fails by −0.032 while our theorem
    still holds. The theorem is strictly stronger and self-contained.
    Check script: `computations/python/gaffke_domination_check.py`.

## Claim 3 — four-coordinate ordered-cap theorem: NO PREDECESSOR FOUND

On 2026-09-08, ten distinct exact-phrase and structural searches covered
Dirichlet averages, B-spline tails, ordered knots, majorization, convex
order, simplex halfspace caps, section centroids, barycentric cuts,
step-vector extrema, Stringer--Gaffke domination, and systematic-PPS
randomization inversion. The 2026-09-17 search added exact identifiers and
structural queries for Gaffke validity, tetrahedral halfspace volumes,
permutation-invariant simplex depth, and weighted Dirichlet tails. No search
proves priority; the conclusion is that no predecessor to the stated
four-coordinate ordered extremum was found in the searched sources.

The classical ingredients and nearest comparisons are distinct:

- Carlson's Dirichlet-average/B-spline theory and Lasserre's Laplace-transform
  derivation give cap and section volume formulas.
- Patáková--Tancer--Wagner, *Discrete & Computational Geometry* 68 (2022),
  Proposition 1.11 (Proposition 11 in the
  [arXiv version](https://arxiv.org/pdf/2003.13536)), characterizes stationary
  halfspace-volume directions by the section centroid. This is a method
  antecedent, not a new stationarity theorem here.
- Paindaveine--Van Bever, *Statistics & Probability Letters* 129 (2017),
  335--339, [Theorem 2](https://arxiv.org/html/1705.04974v1), determines the
  maximal halfspace depth of a class of permutation-invariant simplex
  distributions. For the uniform simplex it covers the equal-barycentric-
  weight special case. It does not give the maximum for arbitrary ordered
  barycentric weights over the matching order cone; its Section 4 leaves
  general off-centroid depth formulas open.
- Webb's 1996 theorem and Myroshnychenko--Tang--Tatarko--Tkocz's 2026
  stability theorem concern central section volume rather than one-sided cap
  volume. Recent barycentric-cut results count suitable hyperplanes rather
  than optimize the ordered cap studied here.
- Diaconis--Perlman, Yu, and Letac--Piccioni provide related weighted-gamma
  crossing and Dirichlet convex-order results, but do not state this ordered
  cap maximum.

The four-coordinate theorem says that for every nondecreasing probability
vector \(c\), the maximum of
\(\Pr\{x\mathbin{\cdot}D>c\mathbin{\cdot}x\}\) over ordered knots \(x\)
is attained at a step vector when \(D\) has four uniform Dirichlet
coordinates. Its new content is the local tetrahedral section inequality,
including the repeated-upper boundary, and its use with the stationary-face
argument. At \(n=3\), this proves pointwise Stringer--Gaffke domination
exactly for
\(0<\alpha\le((19+\sqrt{21})/34)^3=0.3336852118\ldots\).
The endpoint is sharp for this comparison, not asserted to be a Stringer
undercoverage threshold.

### Five- and six-coordinate extensions: CONJECTURAL GLOBAL CONSEQUENCES

The proposed global step-normal theorems in five and six coordinates have
an unresolved proof obligation on repeated-lowest-knot faces. The existing
local computations bound the cap by \(I_{\gamma_0}(1,N-1)\), or by a maximum
that allows all prefixes. The stationary-face argument instead needs the
maximum over *active* prefixes, namely indices \(r\) with
\(x_{r-1}<x_r\). If the lowest knot is repeated, \(r=1\) is not active, and
\(\gamma_0\) is the average of the threshold weights in that block, not
necessarily \(c_0\). Thus the first-coordinate beta bound does not establish
the required active-prefix inequality.

The exact identities, Bernstein multipliers, projective charts, and
six-coordinate family calculations retain their stated local algebraic
conclusions. They do not presently prove the proposed higher-dimensional
global cap theorems. Consequently, the complete \(n=4\) and \(n=5\)
coverage intervals through
\(\alpha_4=0.3214292869970077\ldots\) and
\(\alpha_5=0.3141689898050253\ldots\) are conjectural. The algebraic
Clopper--Pearson weight ranges and terminal-edge comparison obstructions
remain separate results. This identifies a missing implication, not a
counterexample to either global conjecture. Independent direct fixed-level
certificates are not consequences of the unresolved face argument.

### Terminal-edge obstruction

The terminal derivative extends to every sample size. It gives an explicit
critical level \(\alpha_n\) and the limiting constant
\(0.28466813704083846\ldots\), characterized by
\(\alpha(1-2\log\alpha)=1\). No searched source states this
Stringer--Gaffke pointwise-comparison transition. Its proof is the ordinary
one-variable argument in `theory/TERMINAL-EDGE-BIFURCATION.md`; it does not
use either conjectural higher-coordinate global theorem.

## Claim 4 — rigorous fixed-level theorems from n = 5 through n = 7: NO PREDECESSOR FOUND

The 2026 proof of Gaffke's test by Vlassis–Thomas and its bounded-mean
inversion described by Learned-Miller–Thomas and Ming et al. supply a valid
comparison limit. No source found in the review compares that limit
pointwise with Stringer or proves the resulting uniform-simplex cap
inequalities in their stated finite dimensions. No searched source establishes
the corresponding fixed-level \(n=5\), \(n=6\), and \(n=7\) results. The
repository's contribution is the Stringer-specific comparison and rigorous
Bernstein proofs; validity of the Gaffke limit is credited to those sources.
The stronger \(n=3\) continuum is treated in Claim 3. The proposed
\(n=4\) and \(n=5\) continua are conjectural for the reason stated there;
the retained direct certificates establish only their own stated subranges.
No unqualified priority claim is made for these supporting results.

## Claim 5 — all-n Poisson-over-binomial factor transfer: DERIVED FROM A CLASSICAL INEQUALITY

The probability comparison is due to Anderson–Samuels (1967). The supporting
result is to specialize it to the factor roots and use summation by
parts to compare the complete Stringer bounds for every observed sample.
This should be presented as a new audit-specific consequence of a classical
inequality, not as a new probability inequality.

## Claim 6 — direct Poisson simultaneous-band ranges: NO PREDECESSOR FOUND

The proof combines a standard randomized probability integral transform,
the correctly shifted survival-band constraints, and Bolshev's classical
uniform-order-statistic recursion. No reviewed Stringer source applies that
corrected event to the untruncated Poisson audit factors or obtains the
level-dependent ranges \(n\le8\), \(n\le11\), and \(n\le20\). The exact
rational treatment of the Poisson limits and event probabilities is also new
in the reviewed literature. Because the ingredients are classical and the
priority search was not exhaustive, no unqualified first-priority claim is
made for this supporting result.

## Claim 7 — all-n one-cap comparison at confidence at least \(75\%\): NO PREDECESSOR FOUND

The dimension-free cap lemma and the reduction of its
Clopper--Pearson vertices to binomial upper tails are Stringer-specific.
The proof then uses two credited probability results: Anderson--Samuels
(1967) fixed-mean binomial monotonicity and Pinelis's (2021) sharp lower
bound for a binomial exceeding its expectation. No reviewed source combines
these ingredients, proves the resulting factor inequality for all sample
sizes, or derives the zero-uplift region for the Stringer--Gaffke safeguard.
No unqualified priority claim is made for this supporting result.

## Claim 8 — certified two-cap collars: DERIVED EXTENSION, PRIORITY NOT ASSERTED

The two-active-knot spline formula is standard. The contribution is the
ordered-budget prefix-step bound, the one-sign-change argument that reduces a
complete relative-gap collar to its endpoints, and the exact
Clopper--Pearson endpoint certificates through \(n=200\). No reviewed source
was found applying this combination to Stringer, but no unqualified priority
claim is made.

## Claim 9 — search for an earlier correction of Bimpeh's argument

No earlier correction was located in the searched sources. The recorded
citation and the 2002 Bimpeh–Horgan conference abstract do not address the
chapter-5 containment argument. This search result supports only qualified
priority language and is not used to prove the index shift or either
counterexample.

## Claim 10 — first exact explicit counterexample parameters: HOLDS

The only published counterexamples remain Pap–van Zuijlen 1995
(uniform-type taints, existence via recursions, α > 1/2, n > 16); every
later source refers back to them. No prior three-point examples with two
nonzero taint values, no certified rational-arithmetic coverage, and no
violations published in the 30–37% confidence band. Since Pap–van Zuijlen
1996 proved
*asymptotic* conservatism for α ≤ 1/2, certified finite-n violations
approaching the 50% threshold from below sharpen the known picture.

## Claim 11 — systematic-PPS failure and finite-frame inference: QUALIFIED DISTINCTION

Hoogduin, Hall, Tsay, and Pierce (2015), *Auditing: A Journal of Practice &
Theory* 34(4), 85--107, doi:10.2308/ajpt-51081, directly studies whether
systematic MUS evaluated with a simple-random-sample distribution gives
unreliable risk assessments. Its abstract reports simulation evidence of
material risk-assessment error and recommends avoiding that combination. The
companion manuscript therefore does not claim the first warning that
systematic MUS can invalidate a familiar risk calculation.

The mathematical distinction is narrower and stronger: the companion gives
an exact two-item-per-cell family proving binomial Stringer failure for every
\(n\ge2\) and every confidence level, derives conventional-level Poisson
thresholds, proves a necessary diagonal value for every taint-only
universally design-valid rule, constructs a closed-form item-aware overlay,
and specializes exact randomization inversion to a least-favourable-
completion LP that pointwise improves both fallback components. No source
located in the September 8, 2026 search states this failure-and-inference
combination. The 2015 article's complete text was not openly accessible in
this search, and broader systematic-sampling periodicity literature contains
related constructions. Priority language therefore remains qualified; the
theorems rest on the written proofs.

The phase-partition specialization sharpens the distinction. The exact
unequal-phase formula, equal-phase closed form, subset-sum dynamic program,
and direct weak-NP-hardness reduction were not found in the searched audit or
systematic-sampling sources. Fienberg--Neter--Leitch's 1977 multinomial bound
already maximizes an audit estimand over a joint confidence region, and exact
randomization inversion is established methodology in other finite-
population settings. Van Berkel--Vondenhoff (2026) derives exact and
approximate second-order inclusion probabilities for randomized systematic
unequal-probability sampling, but addresses variance estimation rather than
least-favourable taint completion. Kato--Nakagawa (2026) studies exact
sequential finite-population audit boundaries for an attribute-audit design,
not this one-start PPS grid. These sources confirm that the surrounding areas
are active while leaving the fixed-frame phase formula and complexity result
distinct.

Multiple-random-start systematic sampling is longstanding. The companion
does not claim the replication idea. Its narrower result is the exact
all-zero product-design formula and the corresponding unequal-phase dynamic
program for the pre-specified mean-grid ordering.

## July--September 2026 related-work and design-boundary update

- [Vlassis--Thomas, arXiv:2607.08415v1, Theorem 1](https://arxiv.org/html/2607.08415v1),
  proves Gaffke validity for independent nonnegative variables whose individual
  means are at most one. [Ming et al., arXiv:2607.18661v1, Theorem 5.1 and
  Corollary 5.2](https://arxiv.org/html/2607.18661v1), gives the upper endpoint
  as the \(1-\alpha\) quantile of \(D_0+\sum_i t_iD_i\), with the
  finite-sample coverage used here. These exact hypotheses and the endpoint
  representation were checked against the primary texts on 2026-09-17.
  Neither text states the Stringer-specific ordered-cap comparison.
- Martinez-Taboada--Ramdas, arXiv:2608.21694, develops bounded-mean betting
  tests under conditional-mean assumptions and for sampling without
  replacement. That scope is adjacent to, and in one respect broader than,
  the i.i.d. model used here, but it does not establish ordinary Stringer
  coverage or validate systematic PPS sampling through the present proof.
- The repository's SRSWOR result is an elementary conditioning transfer, not
  a claimed new without-replacement concentration inequality. Its contribution
  here is the explicit application to the proved Stringer scopes and to the
  all-sample-size valid alternatives, with the exact no-collision penalty and
  design limitations exposed.
- Berger--Chiodini--Zenga, *Statistical Papers* 62 (2021), 2739--2761,
  explicitly treats systematic unequal-probability MUS without replacement
  and motivates an empirical-likelihood alternative by Stringer's practical
  conservatism. It is now cited to distinguish the paper's safety question
  from efficiency and sampling-design questions.
- Hoogduin--Hall--Tsay--Pierce, *Auditing: A Journal of Practice \& Theory*
  34(4) (2015), 85--107, already supplies simulation evidence that a
  simple-random-sample evaluation can be unreliable after systematic MUS.
  It is now cited as the principal predecessor for the design-risk question.
- The 2025 AICPA audit-sampling guide and PCAOB AS 2315 confirm the continued
  professional relevance of MUS and the broader range of random-based designs.
  They are practice sources, not mathematical support for an i.i.d. theorem.

## Bibliographic note

Łuczak--Mieczkowska--Šileikis, previously recorded as arXiv:1602.03547
(2016), was published in *Statistics & Probability Letters* 129 (2017),
12--16, [doi:10.1016/j.spl.2017.04.024](https://doi.org/10.1016/j.spl.2017.04.024).
The [author's publication list](https://sites.google.com/site/matassileikis/research)
confirms these details. Sources for supporting results need not all appear in
the focused principal article's bibliography.
