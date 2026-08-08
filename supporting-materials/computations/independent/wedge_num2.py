import numpy as np
from itertools import combinations
from wedge_num import params, build_M, pfail, two_atom, three_atom

rng = np.random.default_rng(999)

# 1) print equality-attaining 2-atom configs across w for two alphas
for alpha in [0.5, 0.8]:
    beta, gamma, A, B = params(alpha)
    print(f"alpha={alpha} B={B:.4f} beta={beta:.4f}")
    for w in [0.5*B, B, 0.5*(B+beta), 0.95*beta]:
        f, arg = two_atom(alpha, w, strict=False, ngrid=1201)
        x, y, p = arg
        print(f"  w={w:.4f}: Pfail={f:.8f} atoms x={x:.6f} y={y:.6f} p_x={p:.6f} "
              f"(w/B={w/B:.4f}, (w-B)/A={(w-B)/A:.4f})")

# 2) three-atom exhaustive for a few (alpha, w)
print("\n3-atom exhaustive:")
for alpha in [0.3, 0.5]:
    beta, gamma, A, B = params(alpha)
    for w in [0.7*B, B, 0.5*(B+beta), 0.9*beta]:
        f, arg = three_atom(alpha, w, strict=False, ngrid=25)
        flag = " EXCEEDS" if f > alpha + 1e-7 else ""
        print(f"  alpha={alpha} w={w:.4f}: max={f:.8f} vs alpha={alpha}{flag}")
        if f > alpha + 1e-7:
            print("   ", arg)

# 3) random-support stress: k atoms (k in 2..6), random positions incl. specials,
#    maximize quadratic over the mass polytope via many random vertex-mixture starts
#    + coordinate-pair exchange hill climbing.
def polytope_random_point(xs, w, tries=200):
    n = len(xs)
    for _ in range(tries):
        p = rng.dirichlet(np.ones(n))
        m = p @ xs
        # adjust toward mean w by mixing with a vertex on the right side
        lo = xs.min(); hi = xs.max()
        if not (lo <= w <= hi):
            return None
        # bisection mix with extreme point at argmin/argmax
        j = int(np.argmax(xs)) if m < w else int(np.argmin(xs))
        e = np.zeros(n); e[j] = 1.0
        t = (w - m) / (xs[j] - m) if xs[j] != m else 0.0
        if 0 <= t <= 1:
            return (1 - t) * p + t * e
    return None

def local_max(xs, p, M, w, iters=4000):
    n = len(xs)
    f = pfail(p, M)
    if n < 3:
        return f, p
    for _ in range(iters):
        i, j, k = rng.permutation(n)[:3]
        # move mass among i,j,k preserving sum and mean: direction d with
        # d_i+d_j+d_k=0, d_i x_i + d_j x_j + d_k x_k = 0
        d = np.zeros(n)
        det = (xs[j] - xs[i]) * 1.0
        if abs(xs[k] - xs[j]) < 1e-12 or abs(xs[j] - xs[i]) < 1e-12:
            continue
        # solve d_i + d_j = -d_k ; x_i d_i + x_j d_j = -x_k d_k with d_k = s
        s = (rng.random() - 0.5) * 0.2
        di = (-s * (xs[j] - xs[k])) / (xs[j] - xs[i])
        dj = -s - di
        d[i], d[j], d[k] = di, dj, s
        # max step keeping p >= 0
        neg = d < -1e-15
        tmax = 1.0
        if neg.any():
            tmax = float(np.min(-p[neg] / d[neg]))
        if tmax <= 0:
            continue
        for t in [tmax, 0.5 * tmax, 0.1 * tmax]:
            q = p + t * d
            if (q < -1e-12).any():
                continue
            q = np.maximum(q, 0)
            fq = pfail(q, M)
            if fq > f + 1e-12:
                p, f = q, fq
                break
    return f, p

print("\nrandom-support stress (>= version):")
worst = {}
for alpha in [0.05, 0.3, 0.5, 0.8]:
    beta, gamma, A, B = params(alpha)
    wmax = -1.0; warg = None
    for trial in range(300):
        w = rng.random() * beta * 0.999 + 1e-3
        k = rng.integers(2, 7)
        specials = np.array([0.0, w, min(1, w/beta), min(1, w/B), 1.0])
        pos = np.concatenate([rng.random(k), specials[rng.random(5) < 0.4]])
        xs = np.unique(np.round(pos, 8))
        if len(xs) < 2 or not (xs.min() <= w <= xs.max()):
            continue
        M = build_M(xs, A, B, w, strict=False).astype(float)
        p0 = polytope_random_point(xs, w)
        if p0 is None:
            continue
        f, p = local_max(xs, p0, M, w)
        assert abs(p.sum() - 1) < 1e-6 and abs(p @ xs - w) < 1e-6, (p.sum(), p @ xs, w)
        if f - alpha > wmax:
            wmax = f - alpha; warg = (w, xs, np.round(p, 5), f)
    print(f"  alpha={alpha}: max(Pfail - alpha) over stress = {wmax:.8f}")
    if wmax > 1e-7:
        print("   ", warg)
