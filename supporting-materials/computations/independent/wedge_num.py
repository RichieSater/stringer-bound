"""Independent numerical check of the n=2 wedge inequality.

Claim: alpha in (0,1), beta=sqrt(alpha), gamma=sqrt(1-alpha), A=beta+gamma-1,
B=1-gamma. U1,U2 iid on [0,1], mean w.  P(A min + B max >= w) <= alpha.

We maximize P_fail = sum_{i,j} p_i p_j 1[T(x_i,x_j) >= w] over discrete
distributions with mean w, via:
  (a) exhaustive 2-atom search (analytic mass),
  (b) exhaustive 3-atom search on a structured grid (1-D mass family),
  (c) Frank-Wolfe with exact line search on a fine grid, many restarts.
Both the >= (claim's failure event) and strict > versions are tracked.
"""
import numpy as np
from itertools import combinations

rng = np.random.default_rng(12345)
EPS = 1e-9

def params(alpha):
    beta = np.sqrt(alpha); gamma = np.sqrt(1.0 - alpha)
    return beta, gamma, beta + gamma - 1.0, 1.0 - gamma  # beta,gamma,A,B

def build_M(xs, A, B, w, strict):
    mn = np.minimum.outer(xs, xs); mx = np.maximum.outer(xs, xs)
    T = A * mn + B * mx
    return (T > w + EPS) if strict else (T >= w - EPS)

def pfail(ps, M):
    return float(ps @ M @ ps)

# ---------- (a) two-atom exhaustive ----------
def two_atom(alpha, w, strict, ngrid=2001):
    beta, gamma, A, B = params(alpha)
    # atoms x <= w <= y, mass p at x: p = (y-w)/(y-x)
    special = [0.0, w, min(1.0, w/beta), min(1.0, w/B) if B > 0 else 1.0, 1.0]
    xs = np.unique(np.clip(np.concatenate([np.linspace(0, w, ngrid//2),
                                           np.array([s for s in special if s <= w])]), 0, 1))
    ys = np.unique(np.clip(np.concatenate([np.linspace(w, 1, ngrid//2),
                                           np.array([s for s in special if s >= w])]), 0, 1))
    best = (-1.0, None)
    for x in xs:
        for y in ys:
            if y <= x + 1e-15:
                continue
            p = (y - w) / (y - x)
            if p < -1e-12 or p > 1 + 1e-12:
                continue
            p = min(max(p, 0.0), 1.0)
            Txx = beta * x; Txy = A * x + B * y; Tyy = beta * y
            if strict:
                f = (p*p)*(Txx > w+EPS) + 2*p*(1-p)*(Txy > w+EPS) + ((1-p)**2)*(Tyy > w+EPS)
            else:
                f = (p*p)*(Txx >= w-EPS) + 2*p*(1-p)*(Txy >= w-EPS) + ((1-p)**2)*(Tyy >= w-EPS)
            if f > best[0]:
                best = (f, (x, y, p))
    return best

# ---------- (b) three-atom exhaustive ----------
def three_atom(alpha, w, strict, ngrid=41):
    beta, gamma, A, B = params(alpha)
    special = [0.0, w, w/beta, (w/B if B > 0 else 2.0), 1.0,
               max(0.0, (w - B)/A) if A > 0 else 0.0]
    grid = np.unique(np.clip(np.concatenate([np.linspace(0, 1, ngrid),
                                             np.array([s for s in special if 0 <= s <= 1]),
                                             np.array([s - 1e-6 for s in special if 1e-6 <= s <= 1]),
                                             np.array([s + 1e-6 for s in special if 0 <= s <= 1 - 1e-6])]), 0, 1))
    best = (-1.0, None)
    ts = np.linspace(0, 1, 201)
    for x, y, z in combinations(grid, 3):
        # masses p,q,r >=0, p+q+r=1, px+qy+rz=w  -> 1-D family
        # param by r: p=( (y-w) - r(y-z) )/(y-x), q=1-p-r
        xs3 = np.array([x, y, z])
        M = build_M(xs3, A, B, w, strict).astype(float)
        # r range keeping p,q in [0,1]
        denom = y - x
        if denom < 1e-15:
            continue
        # p(r) = ((y-w) - r*(y-z))/(y-x); q(r) = 1 - p - r
        rmaxcands = []
        feas = []
        for r in ts:
            p = ((y - w) - r * (y - z)) / denom
            q = 1 - p - r
            if p < -1e-12 or q < -1e-12:
                continue
            ps = np.array([max(p, 0), max(q, 0), r]); ps /= ps.sum()
            feas.append(ps)
        if not feas:
            continue
        # objective is quadratic in r; sample + refine at endpoints/critical pt
        vals = [pfail(ps, M) for ps in feas]
        k = int(np.argmax(vals))
        if vals[k] > best[0]:
            best = (vals[k], (x, y, z, tuple(np.round(feas[k], 6))))
    return best

# ---------- (c) Frank-Wolfe on a fine grid ----------
def fw_vertices(xs, w):
    """Vertices of {p>=0, sum p=1, sum p x = w}: <=2-atom distns with mean w."""
    verts = []
    n = len(xs)
    lo = [i for i in range(n) if xs[i] <= w + 1e-12]
    hi = [j for j in range(n) if xs[j] >= w - 1e-12]
    for i in lo:
        for j in hi:
            if abs(xs[j] - xs[i]) < 1e-15:
                if abs(xs[i] - w) < 1e-12:
                    verts.append((i, j, 1.0))
                continue
            p = (xs[j] - w) / (xs[j] - xs[i])
            if -1e-12 <= p <= 1 + 1e-12:
                verts.append((i, j, min(max(p, 0), 1)))
    return verts

def frank_wolfe(alpha, w, strict, n=241, nrestart=40, iters=300):
    beta, gamma, A, B = params(alpha)
    special = [0.0, w, w/beta, (w/B if B > 0 else 2.0), 1.0]
    xs = np.unique(np.clip(np.concatenate([np.linspace(0, 1, n),
                                           np.array([s for s in special if 0 <= s <= 1])]), 0, 1))
    M = build_M(xs, A, B, w, strict).astype(float)
    verts = fw_vertices(xs, w)
    if not verts:
        return (0.0, None)
    best = (-1.0, None)
    N = len(xs)
    for _ in range(nrestart):
        # random feasible start: mix of random vertices
        k = rng.integers(1, 6)
        p = np.zeros(N)
        for _ in range(k):
            i, j, pi = verts[rng.integers(len(verts))]
            lam = rng.random()
            p[i] += lam * pi; p[j] += lam * (1 - pi)
        p /= p.sum()
        # fix mean (mixture of mean-w vertices has mean w already)
        for _ in range(iters):
            g = 2.0 * (M @ p)
            # best vertex for linear obj
            bv, bval = None, -np.inf
            for (i, j, pi) in verts:
                v = pi * g[i] + (1 - pi) * g[j]
                if v > bval:
                    bval, bv = v, (i, j, pi)
            s = np.zeros(N); s[bv[0]] += bv[2]; s[bv[1]] += 1 - bv[2]
            d = s - p
            # f(p+t d) = f(p) + t*(2 p'Md) + t^2 d'Md ; maximize on [0,1]
            a2 = d @ M @ d; a1 = 2 * p @ M @ d
            if a2 <= 0:
                t = 1.0 if a1 + a2 > 0 else (0.0 if a1 <= 0 else min(1.0, -a1/(2*a2)) if a2 < 0 else 1.0)
                # concave: interior max at t=-a1/(2 a2) if in [0,1]
                if a2 < 0:
                    tstar = -a1 / (2 * a2)
                    t = min(max(tstar, 0.0), 1.0)
                else:
                    t = 1.0 if a1 > 0 else 0.0
            else:
                t = 1.0 if a1 + a2 > 0 else 0.0
            if t <= 1e-14:
                break
            p = p + t * d
            p = np.maximum(p, 0); p /= p.sum()
        f = pfail(p, M)
        if f > best[0]:
            supp = [(round(xs[i], 6), round(p[i], 6)) for i in np.nonzero(p > 1e-6)[0]]
            best = (f, supp)
    return best

if __name__ == "__main__":
    for alpha in [0.05, 0.3, 0.5, 0.8]:
        beta, gamma, A, B = params(alpha)
        print(f"\n=== alpha={alpha}  beta={beta:.6f} gamma={gamma:.6f} A={A:.6f} B={B:.6f} ===")
        wlist = sorted(set([0.02, 0.5*B, 0.9*B, B, min(1, 1.05*B), 0.5*(B+beta),
                            0.9*beta, 0.99*beta, beta*0.999, min(beta, 0.999)]))
        wlist = [w for w in wlist if 0 < w < 1]
        for w in wlist:
            f2, w2 = two_atom(alpha, w, strict=False)
            f2s, _ = two_atom(alpha, w, strict=True)
            ffw, sfw = frank_wolfe(alpha, w, strict=False)
            ffws, _ = frank_wolfe(alpha, w, strict=True)
            mx = max(f2, ffw)
            flag = "  <-- EXCEEDS alpha!" if mx > alpha + 1e-7 else ""
            print(f" w={w:.6f}  P>=: 2atom={f2:.8f} FW={ffw:.8f}  (strict: {max(f2s, ffws):.8f})  alpha={alpha}{flag}")
            if mx > alpha + 1e-7:
                print(f"    2atom argmax={w2}  FW support={sfw}")
