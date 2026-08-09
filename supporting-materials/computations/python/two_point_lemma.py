"""Two-point lemma: a single nonzero taint value can never under-cover.

Lemma. Fix ``n``, ``alpha`` and the binomial (Clopper--Pearson) factors
``p_0 <= p_1 <= ... <= p_n``.  Let the taint distribution place probability
``q`` on a single value ``v in (0, 1]`` and ``1 - q`` on 0, so
``theta = q v``.  Then ``P(SB >= theta) >= 1 - alpha``.

Proof.  The number of nonzero taints is ``K ~ Bin(n, q)`` and

    SB = p_0 + v (p_K - p_0) = (1 - v) p_0 + v p_K >= v p_K.

On the event ``{p_K >= q}`` we get ``SB >= v q = theta``.  But
``{p_K >= q}`` is exactly the coverage event of the Clopper--Pearson upper
confidence limit for a binomial proportion, whose coverage is at least
``1 - alpha`` for every ``q`` (classical CP conservatism: ``p_K < q`` iff
``K < k*(q)`` where ``k*(q)`` is the smallest ``j`` with ``p_j >= q``, and
``P(K <= k*(q) - 1) <= alpha`` by the defining equation of ``p_{k*-1}`` and
monotonicity of the binomial tail in ``q``).  QED.

Consequence for the search: the smallest support that can possibly carry a
finite-sample counterexample is ``{v1 > v2 > 0}`` plus an atom at 0, which
is why ``search_two_value.py`` starts there.

Poisson factors.  The same argument gives conservatism for the Poisson
(AICPA) factors whenever ``p_j^pois >= p_j^binom`` for all ``j``, since then
``SB_pois >= SB_binom`` pointwise on every sample.  That domination is not
proved here in general; this module checks it at high precision for every
``n`` in the requested range, and separately brute-forces the lemma's
conclusion on a grid as an independent test of the coverage machinery.

Usage:
    python3 two_point_lemma.py --alpha 0.05 --n-max 100
"""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction

from mpmath import mp, mpf

from coverage import coverage_exact
from stringer import factors


def check_poisson_dominates(n: int, alpha: str, dps: int = 50) -> bool:
    """Numerically check separated factor brackets for all ``0 <= j <= n``.

    This is a high-precision check, not a directed-rounding or interval-
    arithmetic proof of the tail evaluations.
    """
    fb = factors(n, alpha, "binomial", dps)
    fp = factors(n, alpha, "poisson", dps)
    with mp.workdps(dps):
        for j in range(n + 1):
            pb, wb = fb[j]
            pp, wp = fp[j]
            if not (pp - wp > pb + wb):
                return False
    return True


def brute_force_two_point(n: int, alpha: str, v_steps: int = 20,
                          q_steps: int = 20, dps: int = 50):
    """Exact coverage over a (v, q) grid for the single-value support.

    Returns the minimum exact coverage found.  By the lemma this must be
    >= 1 - alpha; a smaller value would indicate a bug in the machinery.
    """
    worst = None
    for iv in range(1, v_steps + 1):
        v = Fraction(iv, v_steps)
        for iq in range(1, q_steps):
            q = Fraction(iq, q_steps)
            cov, _theta, _margin = coverage_exact([v], [q], n, alpha,
                                                  "binomial")
            if worst is None or cov < worst[0]:
                worst = (cov, v, q)
    return worst


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", default="0.05")
    ap.add_argument("--n-max", type=int, default=100)
    ap.add_argument("--brute-n", nargs="*", type=int, default=[5, 10, 20],
                    help="sample sizes for the exact brute-force cross-check")
    args = ap.parse_args(argv)

    nominal = 1 - Fraction(args.alpha)

    print("Poisson-dominates-binomial high-precision factor check:")
    bad = [n for n in range(1, args.n_max + 1)
           if not check_poisson_dominates(n, args.alpha)]
    if bad:
        print("  DOMINATION FAILS for n in %r -- the Poisson lemma does NOT "
              "follow for these n" % bad)
    else:
        print("  p_j^poisson > p_j^binomial for all j <= n, all n <= %d: OK"
              % args.n_max)

    print("Exact brute-force of the lemma conclusion (binomial factors):")
    ok = True
    for n in args.brute_n:
        cov, v, q = brute_force_two_point(n, args.alpha)
        status = "OK" if cov >= nominal else "VIOLATION (machinery bug?)"
        if cov < nominal:
            ok = False
        print("  n=%-4d min exact coverage on grid = %s ~ %.10f at v=%s q=%s"
              "  [%s]" % (n, cov, float(cov), v, q, status))
    if not ok or bad:
        return 1
    print("two_point_lemma: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
