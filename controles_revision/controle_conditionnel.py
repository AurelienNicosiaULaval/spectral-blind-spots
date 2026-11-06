"""Independent simulation of the canonical projected sign-flip diagnostic under population-blind
phase alternatives. Written from the manuscript (Section 3), without the authors' R code.
Usage: python3 -I controle_conditionnel.py OUTPUT_DIR"""
import sys, os, json, math
import numpy as np
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
TWO_PI = 2 * np.pi
ORIG = [(0, 0), (0, np.pi / 2), (np.pi / 2, 0), (np.pi / 2, np.pi / 2)]
K = np.arange(0, 9); L = np.arange(1, 9)

def draw_y(x, order, eta, rng):
    y = np.empty_like(x); todo = np.arange(len(x))
    while todo.size:
        prop = rng.uniform(0, TWO_PI, todo.size)
        acc = rng.uniform(size=todo.size) < (1 + eta * np.sin(order * (x[todo] - prop))) / (1 + eta)
        y[todo[acc]] = prop[acc]; todo = todo[~acc]
    return y

def contributions(x, y, z):
    rows = []
    for a, b in ORIG:
        X = np.cos(np.outer(K, x - a))                         # 9 x n
        D = np.cos(np.outer(L, z - b)) - np.cos(np.outer(L, y - b))  # 8 x n
        rows.append((X[:, None, :] * D[None, :, :]).reshape(-1, len(x)))
    return np.vstack(rows)                                    # 288 x n

def design_means(x, order, eta):
    """Exact conditional means E[c_i(t) | x_i] for the uniform candidate."""
    rows = []
    for a, b in ORIG:
        X = np.cos(np.outer(K, x - a))
        D = np.zeros((len(L), len(x)))
        if order in L:
            D[list(L).index(order)] = -(eta / 2) * np.sin(order * (x - b))
        rows.append((X[:, None, :] * D[None, :, :]).reshape(-1, len(x)))
    return np.vstack(rows)

def sign_test(C, S):
    den = np.sqrt((C ** 2).sum(1)); act = den > den.max() * 1e-12
    Cs = C[act] / den[act, None]
    obs = np.abs(Cs.sum(1)).max(); ran = np.abs(Cs @ S).max(0)
    return (1 + np.sum(ran >= obs - 1e-12 * max(1, obs))) / (S.shape[1] + 1)

def run(order, eta, n, R, seed, null=False):
    rng = np.random.default_rng(seed); rej = []; zeta = []
    for r in range(R):
        x = rng.uniform(0, TWO_PI, n)
        y = draw_y(x, order, eta, rng)
        z = draw_y(x, order, eta, rng) if null else rng.uniform(0, TWO_PI, n)
        S = rng.choice([-1.0, 1.0], size=(n, 199))
        C = contributions(x, y, z)
        rej.append(sign_test(C, S) <= 0.05)
        if not null:
            M = design_means(x, order, eta)
            zeta.append(np.max(np.abs(M.sum(1)) / np.sqrt((C ** 2).sum(1))))
    rej = np.array(rej)
    k = int(rej.sum()); p = k / R; zc = 1.959963984540054
    c = (p + zc**2/(2*R)) / (1 + zc**2/R); w = zc*math.sqrt(p*(1-p)/R + zc**2/(4*R*R)) / (1 + zc**2/R)
    res = dict(order=order, eta=eta, n=n, R=R, rejections=k, rate=p, wilson=[c - w, c + w])
    if zeta:
        zeta = np.array(zeta); q = np.quantile(zeta, [0.25, 0.5, 0.75])
        res["max_standardized_design_mean_quartiles"] = [float(v) for v in q]
        hi = zeta > np.median(zeta)
        res["rate_if_design_mean_above_median"] = float(rej[hi].mean())
        res["rate_if_design_mean_below_median"] = float(rej[~hi].mean())
    return res

def fixed_designs(order, eta, n, designs, R, seed):
    rng = np.random.default_rng(seed); rates = []; zetas = []
    for d in range(designs):
        x = rng.uniform(0, TWO_PI, n)
        M = design_means(x, order, eta).sum(1)
        rej = 0
        for r in range(R):
            y = draw_y(x, order, eta, rng); z = rng.uniform(0, TWO_PI, n)
            S = rng.choice([-1.0, 1.0], size=(n, 199))
            rej += sign_test(contributions(x, y, z), S) <= 0.05
        rates.append(rej / R)
        # sd of c_i(t): approximate population value for standardization
        zetas.append(float(np.max(np.abs(M)) / math.sqrt(n * 0.5)))
    return dict(order=order, eta=eta, n=n, designs=designs, replicates_per_design=R,
                rates=[float(v) for v in rates], approx_max_noncentrality=zetas)

res = {"note": "Independent Python implementation of the canonical diagnostic (four sine/cosine views, moments k=0..8, l=1..8, 199 sign vectors).",
       "unconditional": [], "fixed_design": []}
for (order, eta, null) in [(1, 0.0, True), (4, 0.8, False), (2, 0.8, False), (2, 0.4, False), (6, 0.8, False)]:
    for n in (240, 720):
        r = run(order, eta, n, 2000, 1000 * order + n + int(10 * eta), null=null)
        r["null"] = null; res["unconditional"].append(r); print(r, flush=True)
res["fixed_design"].append(fixed_designs(4, 0.8, 240, 12, 400, 777))
print(res["fixed_design"][-1], flush=True)
json.dump(res, open(os.path.join(out, "controle_conditionnel.json"), "w"), indent=1)
