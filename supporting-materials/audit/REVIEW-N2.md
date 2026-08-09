# Review record: the n = 2 theorem

This file records the adversarial reviews of
[`theory/N2-PROOF.md`](../theory/N2-PROOF.md). The reviews are aids to
checking; they are not substitutes for independent human peer review, and
the written proof remains the authority.

## Initial adversarial passes

Two separate AI-assisted passes attacked the factor formulas, reduction,
rectangle lemma, boundary limits, potential inequality, change of variables,
strictness cases, and sharpness construction. They also ran structured and
random searches over atomic distributions. Those passes found no violation
of the conservatism inequality. Their useful corrections included explicit
treatment of \(w=0\), the convention \(h_+(1)=0\), endpoint checks in both
change-of-variable regimes, and a more precise statement of the two
sharpness families.

One pass also developed an alternative derivation. Because that derivation
was AI-assisted rather than independently refereed by another researcher, it
is treated as corroboration only and is not used in the manuscript.

## Subsequent mathematical audit at commit `3c6c3a2`

A later line-by-line audit identified a genuine gap in the initial proof of
the theorem's **non-attainment** clause. The weak equality argument had
claimed that the potential satisfies equality only at \(g=\beta\). In fact,

\[
\widetilde k(B)=\widetilde k(\beta)=\beta.
\]

The strict conservatism proof excludes \(g=B\) through its pointwise bound
\(G>B\), so that proof was unaffected. The weak attainment argument gives
only \(G\ge B\), however, and therefore had to retain both equality values.

The proof has been repaired as follows. Equality in the budget chain and
monotonicity reduce the survival function to

\[
G(x)=
\begin{cases}
\beta,&0\le x<c,\\
B,&c\le x<d,\\
0,&d\le x\le1,
\end{cases}
\qquad Ac+Bd=w.
\]

Every such equality case assigns positive probability to
\(Am+BM=w\): the mixed pair \((c,d)\) when \(0<c<d\), the mixed pair
\((0,d)\) when \(c=0<d\), or the pair \((u^*,u^*)\) when \(c=d=u^*\).
Attainment would require this boundary probability to be zero, giving the
missing contradiction.

The same audit identified and prompted these additional corrections:

- deletion of an invalid scaling parenthetical that had claimed equivalence
  between the weak and strict wedge forms;
- recognition that \(A=B\) at the interior value
  \(\alpha=16/25\), rather than \(A\ne B\) for every \(\alpha\);
- qualification of numerical Poisson-factor domination claims; and
- replacement of “decreasing weights” by “positive increments” in the
  general-\(n\) discussion.

## Current verification status

The audit reported that the following parts passed line-by-line checking:
the explicit factors, wedge reduction, boundary regimes, rectangle and
survival bounds, boundary-pair algebra, potential formulas and concavity,
integral contradiction, sharpness construction, and the proposition for one
nonzero taint value. The repaired non-attainment argument is now printed in
both the manuscript and `N2-PROOF.md`.

`n2_proof_check.py` symbolically checks the algebraic identities and stress
tests the theorem over approximately 200,000 atomic distributions. Those
computations can expose implementation or algebra errors; they do not replace
the proof.
