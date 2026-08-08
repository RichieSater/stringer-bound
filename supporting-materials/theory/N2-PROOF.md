# Theorem: the Stringer bound is conservative at n = 2

**Theorem.** Let \(\alpha \in (0,1)\) and let \(p_0 = 1-\sqrt{\alpha}\),
\(p_1 = \sqrt{1-\alpha}\), \(p_2 = 1\) be the binomial (Clopper–Pearson)
Stringer factors for \(n = 2\). For every distribution \(F\) on \([0,1]\)
and \(T_1, T_2\) i.i.d. \(\sim F\) with \(\mu = E[T]\),

\[ P(\mathrm{SB} \ge \mu) \;\ge\; 1-\alpha . \]

The bound is sharp: the supremum of \(P(\mathrm{SB} < \mu)\) over all
\(F\) equals \(\alpha\), approached but not attained.

**Remarks.** (1) This is, to our knowledge, the first finite-sample proof
of Stringer conservatism at any \(n \ge 2\) free of hypotheses on \(F\)
(Bickel 1992 proved \(n = 1\); his \(n \ge 2\) statement gives only
\((1-\alpha)^{n+1}\) "under certain conditions"). (2) The theorem holds
for **all** \(\alpha \in (0,1)\), including \(\alpha > 1/2\) where the
bound is asymptotically anti-conservative — so the known violations
(Pap–van Zuijlen 1995; exact certificates in this repository at
\(n \ge 50\)) are strictly a large-\(n\) phenomenon. (3) Poisson-factor
version follows a fortiori wherever the Poisson factors dominate the
binomial ones (checked rigorously per \(n\) in `two_point_lemma.py`).

Throughout write \(\beta = \sqrt{\alpha}\), \(\gamma = \sqrt{1-\alpha}\)
(so \(\beta^2 + \gamma^2 = 1\)), and

\[ A = p_1 - p_0 = \beta + \gamma - 1 \;>\; 0, \qquad
   B = 1 - p_1 = 1 - \gamma \;>\; 0, \qquad A + B = \beta . \]

(\(A > 0\) because \((\beta+\gamma)^2 = 1 + 2\sqrt{\alpha(1-\alpha)} > 1\).)

## Reduction

Set \(U_i = 1 - T_i \in [0,1]\), \(w = 1 - \mu = E[U]\),
\(m = \min(U_1,U_2)\), \(M = \max(U_1,U_2)\). Substituting into the
Stringer bound and telescoping gives the identity (machine-checked in
`n2_proof_check.py`)

\[ \mathrm{SB} \;=\; 1 - \big(A\,m + B\,M\big), \]

so \(\mathrm{SB} \ge \mu \iff A m + B M \le w\). The theorem is
implied by: **for \(U_1,U_2\) i.i.d. on \([0,1]\) with mean
\(w > 0\),** (the \(w = 0\) case is handled separately below — note
(&#42;) itself is false at \(w = 0\), where the strict event is empty)

\[ P\big(A m + B M < w\big) \;\ge\; \gamma^2 . \tag{&#42;} \]

(Indeed \(P(\mathrm{SB} \ge \mu) = P(Am+BM \le w) \ge P(Am+BM<w)\).
The strict-inequality form (&#42;) is in fact equivalent to the theorem —
apply the weak form to \((1-\varepsilon)U\) and let
\(\varepsilon \downarrow 0\) — but only the stated implication is used.)

Degenerate regimes. If \(w = 0\) then \(U \equiv 0\) a.s. and
\(\mathrm{SB} = 1 \ge \mu\) surely. If \(w > \beta\) the event
\(\{Am+BM \ge w\}\) is empty (its maximum over the square is
\(A+B = \beta\)), and if \(w = \beta\) it is \(\{U_1 = U_2 = 1\}\), of
probability \(P(U=1)^2 \le w^2 = \beta^2\) by Markov; either way
\(P(Am+BM \ge w) \le \alpha\), which is (&#42;) and more. **Assume from
now on \(0 < w < \beta\).**

## Setup

Let \(F(x) = P(U \le x)\), \(G = 1 - F\), and recall
\(w = E[U] = \int_0^1 G(t)\,dt\). Define

\[ u^* = \frac{w}{\beta} \in (0,1), \qquad
   a_0 = \max\Big(0, \frac{w-B}{A}\Big), \qquad
   T = \min\Big(\frac{w}{B}, 1\Big). \]

Elementary checks (using \(0 < w < \beta\) and \(A, B > 0\)):
\(a_0 < u^* < T\).

**Rectangle lemma.** If \(0 \le a \le b \le 1\) and \(Aa + Bb < w\), then

\[ P(Am+BM < w) \;\ge\; 2F(a)F(b) - F(a)^2 . \]

*Proof.* On the event \(\{U_1 \le a, U_2 \le b\} \cup \{U_1 \le b, U_2
\le a\}\) we have \(m \le a\) and \(M \le b\), hence \(Am + BM \le Aa +
Bb < w\) (using \(A, B \ge 0\)). The two events have union probability
\(2F(a)F(b) - F(a)^2\) by inclusion–exclusion. ∎

Assume now, for contradiction, that (&#42;) fails:

\[ \textbf{(H)} \qquad 2F(a)F(b) - F(a)^2 \;<\; \gamma^2
   \quad\text{for all } 0 \le a \le b \le 1 \text{ with } Aa+Bb < w. \]

## Three families of pointwise bounds under (H)

**(i) Squares.** For \(t < u^*\): \(a = b = t\) is admissible
(\(\beta t < w\)), and the lemma value is \(F(t)^2\). So \(F(t)^2 <
\gamma^2\), i.e.

\[ G(t) > 1 - \gamma = B \qquad \text{for all } t \in [0, u^*). \]

**(ii) Full-height rectangles** (nonvacuous iff \(a_0 > 0\), i.e.
\(w > B\)). For \(a < a_0\): \((a, 1)\) is admissible (\(Aa + B < w\)),
value \(2F(a) - F(a)^2 = 1 - G(a)^2\). So \(1 - G(a)^2 < \gamma^2\), i.e.

\[ G(a) > \beta \qquad \text{for all } a \in [0, a_0). \]

**(iii) Boundary pairs.** Define \(h(g) = \dfrac{\beta^2 - g^2}{2(1-g)}\)
for \(g \in [0, 1)\). For \(b \in (u^*, T)\) let
\(\bar a(b) = (w - Bb)/A \in (a_0, u^*)\). For every \(a' < \bar a(b)\)
the pair \((a', b)\) is admissible; letting \(a' \uparrow \bar a(b)\) in
(H) and writing \(g = G(\bar a(b)^-)\):

- if \(g = 1\) the conclusion below is vacuous and not needed;
- otherwise \(F(\bar a^-) = 1-g > 0\) and (H) gives
  \(F(b) \le \dfrac{\gamma^2 + (1-g)^2}{2(1-g)}\), i.e.

\[ G(b) \;\ge\; h\big(G(\bar a(b)^-)\big) \qquad
   \text{for all } b \in (u^*, T). \]

## The potential function

Let \(h_+ = \max(h, 0)\), extended to \(g = 1\) by \(h_+(1) := 0\)
(note \(h(g) \to -\infty\) as \(g \uparrow 1\), so \(h_+\) is continuous
on \([0,1]\); this convention covers the \(F(\bar a^-) = 0\) case of
(iii), where the bound \(G(b) \ge 0 = h_+(1)\) is what the budget below
actually uses). Define

\[ \tilde k(g) \;=\; g + \frac{A}{B}\, h_+(g), \qquad g \in (B, 1]. \]

**Claim: \(\tilde k(g) \ge \beta\) on \((B, 1]\), with equality only at
\(g = \beta\).** Verified facts (symbolic, `n2_proof_check.py`):

- \(\beta^2 - B^2 = 2\gamma(1-\gamma)\), whence \(h(B) = B\) and so
  \(\tilde k(B) = B + \tfrac{A}{B}B = A + B = \beta\);
- \(h(\beta) = 0\), so \(\tilde k(\beta) = \beta\);
- \(h''(g) = -\gamma^2/(1-g)^3 < 0\) on \([0,1)\): \(h\) is strictly
  concave; since \(h(B) = B > 0\) and \(h(\beta) = 0\) and \(h\) is
  decreasing past its critical point, \(h \ge 0\) on \([B, \beta]\) and
  \(h < 0\) on \((\beta, 1)\).

On \([B, \beta]\): \(\tilde k = g + \tfrac{A}{B}h(g)\) is strictly
concave with equal endpoint values \(\beta\), hence \(\ge \beta\) there
and \(> \beta\) strictly on the open interval \((B, \beta)\). On
\((\beta, 1]\): \(h_+ = 0\), so \(\tilde k(g) = g > \beta\). ∎

## Budget contradiction

Using \(w = \int_0^1 G\), split at \(a_0\), \(u^*\), \(T\) and discard
\([T, 1]\):

\[ w \;\ge\; \int_0^{a_0} G \;+\; \int_{a_0}^{u^*} G \;+\;
   \int_{u^*}^{T} G . \]

For the tail, substitute \(t = b(a) = (w - Aa)/B\), which maps
\([a_0, u^*]\) decreasingly onto \([u^*, T]\) with \(|dt| =
(A/B)\,da\), and apply (iii) (noting \(G(a^-) = G(a)\) for all but
countably many \(a\), so the integrals agree):

\[ \int_{u^*}^{T} G(t)\,dt \;=\; \frac{A}{B}\int_{a_0}^{u^*}
   G(b(a))\,da \;\ge\; \frac{A}{B}\int_{a_0}^{u^*} h_+\big(G(a)\big)\,da . \]

By (i), \(G > B\) on \([a_0, u^*)\), so the potential bound applies
pointwise:

\[ \int_{a_0}^{u^*}\Big[G(a) + \tfrac{A}{B}h_+(G(a))\Big] da
   \;=\; \int_{a_0}^{u^*} \tilde k(G(a))\,da \;\ge\; \beta\,(u^*-a_0). \]

By (ii), \(G > \beta\) pointwise on \([0, a_0)\), so \(\int_0^{a_0} G
\ge \beta a_0\), **strictly** if \(a_0 > 0\). Combining:

\[ w \;\ge\; \beta a_0 + \beta (u^* - a_0) \;=\; \beta u^* \;=\; w, \]

with equality forced throughout. It remains to exhibit one strict
inequality in every case:

- **If \(a_0 > 0\)** (\(w > B\)): \(\int_0^{a_0} G > \beta a_0\) is
  strict (a measurable function exceeding \(\beta\) pointwise on a set
  of positive measure integrates to strictly more). Contradiction.
- **If \(a_0 = 0\)** (\(w \le B\)) **and \(G \ne \beta\) on a
  positive-measure subset of \((0, u^*)\)**: there \(\tilde k(G) >
  \beta\) strictly (Claim), so \(\int_0^{u^*}\tilde k(G) > \beta u^*\).
  Contradiction.
- **If \(a_0 = 0\) and \(G = \beta\) a.e. on \((0, u^*)\)**: for each
  \(b \in (u^*, T)\) choose \(a' < \bar a(b)\) with \(G(a') = \beta\)
  (possible a.e.); the strict inequality (H) at the admissible corner
  \((a', b)\) reads \((1-\beta)\big(1 + \beta - 2G(b)\big) <
  \gamma^2 = (1-\beta)(1+\beta)\), which forces \(G(b) > 0\) for
  **every** \(b \in (u^*, T)\). Hence \(\int_{u^*}^{T} G > 0\) and
  \(w \ge \int_0^{u^*}G + \int_{u^*}^T G > \beta u^* = w\).
  Contradiction.

All cases contradict \(w = w\). So (H) is impossible, proving (&#42;)
and the theorem. ∎

## Sharpness

Two families approach equality (exact case analysis, machine-checked):

- *Diagonal*: mass \(q\) at \(u\) and \(1-q\) at \(0\) (taints
  \(\{1, 1-u\}\), mean \(w = qu\)): the failure event is
  \(\{U_1 = U_2 = u\}\) precisely when \(\beta u > qu\), i.e.
  \(q < \beta\); failure probability \(q^2 \uparrow \beta^2 = \alpha\)
  as \(q \uparrow \beta\).
- *Wing*: mass \(q\) at \(u\) and \(1-q\) at \(0\) with \(q < B\): the
  pair \(\{0, u\}\) fails precisely when \(Bu > qu\), i.e. \(q < B\),
  giving failure probability \(1-(1-q)^2 \uparrow 1-\gamma^2 = \alpha\)
  as \(q \uparrow B\).

A one-parameter three-atom family interpolates between them with
failure probability \(< \alpha\) throughout and \(= \alpha\) in both
limits.

**Non-attainment.** No \(F\) attains \(P(\mathrm{SB} < \mu) = \alpha\).
Suppose some \(F\) did; then \(P(Am+BM > w) = \alpha\), so
\(P(Am+BM < w) = \gamma^2 - P(Am+BM = w) \le \gamma^2\), and the weak
form of (H) holds (every rectangle bound is \(\le \gamma^2\)). Rerunning
the budget chain with weak inequalities forces equality in every link:
equality in \(\int_0^{a_0} G = \beta a_0\) (where now \(G \ge \beta\)
pointwise) and in \(\int_{a_0}^{u^*}\tilde k(G) = \beta(u^*-a_0)\)
forces \(G = \beta\) a.e. on \((0, u^*)\) regardless of the value of
\(a_0\) (since \(\tilde k > \beta\) strictly off \(g = \beta\), which
survives weakening); equality in the tail forces
\(G = h_+(\beta) = 0\) a.e. on \((u^*, T)\) and \(\int_T^1 G = 0\). A
right-continuous monotone \(G\) equal to \(\beta\) a.e. on \((0,u^*)\)
and \(0\) a.e. beyond is exactly \(F = (1-\beta)\delta_0 +
\beta\delta_{u^*}\). But for that \(F\) the largest value of
\(Am + BM\) is \(\beta u^* = w\), attained by the pair \((u^*, u^*)\),
so \(P(Am+BM > w) = 0 \ne \alpha\). Contradiction. ∎

## Verification

`n2_proof_check.py` re-checks: the reduction identity, every symbolic
identity and the concavity used above (sympy), \(\tilde k \ge \beta\)
numerically across \(\alpha\), and the theorem statement itself against
~60,000 atomic distributions per confidence level (max observed failure
probability minus \(\alpha\): exactly \(0\) at the boundary
configurations, negative elsewhere; never positive).
