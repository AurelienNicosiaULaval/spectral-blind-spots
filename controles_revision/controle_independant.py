"""Independent targeted checks for 'Spectral Blind Spots in Circular Imputation Diagnostics'.
Written from the manuscript statements only; it does not import or call the authors' R code.
Usage: python3 -I controle_independant.py OUTPUT_DIR
"""
import sys, json, math, os
import numpy as np
from scipy import special, integrate, stats

out = sys.argv[1]; os.makedirs(out, exist_ok=True)
rep = {}
TWO_PI = 2 * np.pi

# ---------- 1. Null space of the four views at common origins (grid-based, no block algebra) ----------
def basis(N):
    """Real tensor trigonometric basis on T^2, orders 0..N, without the constant."""
    fs = []
    for k in range(N + 1):
        for l in range(N + 1):
            cx = [np.cos] if k == 0 else [np.cos, np.sin]
            cy = [np.cos] if l == 0 else [np.cos, np.sin]
            for f in cx:
                for g in cy:
                    if k == 0 and l == 0:
                        continue
                    fs.append((k, l, f, g))
    return fs

def view_gram(N, origins, G=48, crossed=None):
    """Gram matrix of the linear map coefficients -> concatenated double symmetrizations on a grid.
    common origins: views (X-a-j*pi/2, Theta-a-m*pi/2), j,m in {0,1} (cosines of these are the sine/cosine views).
    crossed: list of (a,b) pairs, cosine views only. Accumulated view by view to save memory."""
    u = TWO_PI * (np.arange(G) + 0.37) / G
    U, V = np.meshgrid(u, u, indexing="ij")
    fs = basis(N)
    pairs = []
    if crossed is None:
        for a in origins:
            for j in (0, 1):
                for m in (0, 1):
                    pairs.append((a + j * np.pi / 2, a + m * np.pi / 2))
    else:
        pairs = crossed
    Gm = np.zeros((len(fs), len(fs)))
    for (a, b) in pairs:
        cols = []
        for (k, l, f, g) in fs:
            xs = [U, 2 * a - U]; ts = [V, 2 * b - V]
            s = sum(f(k * x) * g(l * t) for x in xs for t in ts) / 4.0
            cols.append(s.ravel())
        Av = np.column_stack(cols)
        Gm += Av.T @ Av
    return Gm, fs

def null_info(Gm, rel=1e-12):
    ev, vec = np.linalg.eigh(Gm)
    keep = ev < rel * ev.max()
    return int(keep.sum()), vec[:, keep], float(np.sqrt(max(ev[~keep].min(), 0) / ev.max()))

N = 10
origins_all = TWO_PI * np.arange(41) / 41
Gm, fs = view_gram(N, origins_all)
dim_null, Nul, rel_gap = null_info(Gm)
expected = []
for k in range(2, N + 1, 2):
    v = np.zeros(len(fs))
    for i, (kk, ll, f, g) in enumerate(fs):
        if kk == k and ll == k and f is np.sin and g is np.cos: v[i] = 1
        if kk == k and ll == k and f is np.cos and g is np.sin: v[i] = -1
    expected.append(v / np.linalg.norm(v))
E = np.array(expected)
cosines = np.linalg.svd(E @ Nul, compute_uv=False) if dim_null else np.array([0.0])
rep["nullspace_common_all_origins"] = dict(order=N, n_basis=len(fs), n_origins=len(origins_all),
    null_dimension=dim_null, expected_dimension=len(expected),
    min_cos_principal_angle=float(cosines.min()), relative_smallest_nonnull_singular_value=rel_gap)

# finite set of four common origins used in the simulation, orders <= 8
G4, fs8 = view_gram(8, [0, np.pi / 8, np.pi / 3, np.pi / 4])
d4, N4, g4 = null_info(G4)
extra = [ (fs8[i][0], fs8[i][1]) for i in np.unique(np.nonzero(np.abs(N4) > 1e-6)[0]) ]
rep["nullspace_four_common_origins_order8"] = dict(n_basis=len(fs8), null_dimension=d4,
    expected_if_no_extra_blind_spots=4, blocks_involved=sorted(set(extra)), relative_smallest_nonnull_singular_value=g4)
# crossed pi/7, orders <=6 should be injective; order 7 adds resonant directions
for name, ang, M in [("pi7_order6", np.pi / 7, 6), ("pi7_order7", np.pi / 7, 7), ("pi4_order3", np.pi / 4, 3), ("pi4_order4", np.pi / 4, 4)]:
    Gc, fsc = view_gram(M, None, crossed=[(0, 0), (0, ang), (ang, 0), (ang, ang)])
    dc, Nc, gc = null_info(Gc)
    blocks = sorted(set((fsc[i][0], fsc[i][1]) for i in np.unique(np.nonzero(np.abs(Nc) > 1e-6)[0]))) if dc else []
    rep["crossed_" + name] = dict(null_dimension=dc, blocks_involved=blocks, relative_smallest_nonnull_singular_value=gc)

# ---------- 2. Maximin angle by brute force ----------
mm = {}
grid = np.linspace(1e-6, np.pi - 1e-6, 400001)
for M in range(1, 11):
    k = np.arange(1, M + 1)[:, None]
    val = np.sqrt(1 - np.abs(np.cos(k * grid))).min(axis=0)
    best = val.max(); arg = grid[np.isclose(val, best, atol=1e-6)]
    # cluster maximizers
    clusters = []
    for a in arg:
        if not clusters or a - clusters[-1][-1] > 1e-3: clusters.append([a])
        else: clusters[-1].append(a)
    centers = [float(np.mean(c) / np.pi * (M + 1)) for c in clusters]
    mm[M] = dict(best=float(best), formula=float(np.sqrt(1 - np.cos(np.pi / (M + 1)))),
                 maximizers_in_units_of_pi_over_M_plus_1=[round(c, 3) for c in centers])
rep["maximin"] = mm

# ---------- 3. Phase contrast formula ----------
def proj_contrast(k, eta, a, b):
    f = lambda t, x: np.cos(k * (x - a)) * np.cos(k * (t - b)) * eta * np.sin(k * (x - t)) / (4 * np.pi ** 2)
    v = integrate.dblquad(f, 0, TWO_PI, 0, TWO_PI, epsabs=1e-11)[0]
    return -v  # E_Q - E_P with Q uniform
errs = []
rng = np.random.default_rng(20261010)
for _ in range(6):
    k = int(rng.integers(1, 8)); a, b = rng.uniform(0, TWO_PI, 2)
    errs.append(abs(proj_contrast(k, 0.8, a, b) - (-0.8 / 4 * np.sin(k * (a - b)))))
rep["phase_contrast_formula_max_error"] = float(max(errs))

# ---------- 4. Arc risk and circular CRPS identities ----------
def arc(a, b): return np.abs(np.angle(np.exp(1j * (a - b))))
G = 4096; t = TWO_PI * np.arange(G) / G; w = TWO_PI / G
def vm_mix(mus, kap, wts):
    d = sum(wt * np.exp(kp * np.cos(t - m)) / (TWO_PI * special.i0(kp)) for m, kp, wt in zip(mus, kap, wts))
    return d / (d.sum() * w)
P = vm_mix([0.3, 2.9], [3.0, 8.0], [0.6, 0.4]); Q = vm_mix([1.1], [1.5], [1.0])
D = arc(t[:, None], t[None, :])
def mom(p, k): return np.sum(p * np.exp(1j * k * t)) * w
z0 = 0.7
lhs = np.sum(P * arc(t, z0)) * w
rhs = np.pi / 2 - 4 / np.pi * sum((mom(P, k) * np.exp(-1j * k * z0)).real / k ** 2 for k in range(1, 4000, 2))
def crps_expect(Qd, Pd):  # E_{Y~P} S(Q,Y)
    return (Pd @ D @ Qd) * w * w - 0.5 * (Qd @ D @ Qd) * w * w
excess = crps_expect(Q, P) - crps_expect(P, P)
formula = 2 / np.pi * sum(abs(mom(P, k) - mom(Q, k)) ** 2 / k ** 2 for k in range(1, 4000, 2))
rep["arc_identity_abs_error"] = float(abs(lhs - rhs))
rep["crps_identity"] = dict(excess=float(excess), formula=float(formula), abs_error=float(abs(excess - formula)))
# antipodal law: arc risk pi/2 for every z
Pa = (1 + 0.8 * np.sin(2 * (0.23 - t))) / TWO_PI
rep["antipodal_arc_risk_range"] = [float(min(np.sum(Pa * arc(t, z)) * w for z in t[::256])), float(max(np.sum(Pa * arc(t, z)) * w for z in t[::256]))]

# ---------- 5. Kernel Fourier coefficients ----------
c = 1 / 0.5 ** 2
kern = np.exp(-(1 - np.cos(t)) / 0.5 ** 2)
fft = np.fft.fft(kern).real / G
coef = np.exp(-c) * special.iv(np.arange(0, 13), c)
rep["kernel_coefficients"] = dict(max_abs_error=float(np.max(np.abs(fft[:13] - coef))),
                                  coefficients_0_to_8=[float(x) for x in coef[:9]],
                                  ratio_order4_to_order1=float(coef[4] / coef[1]),
                                  ratio_order6_to_order1=float(coef[6] / coef[1]))

# ---------- 6. Population distances ----------
bump = lambda x, c0: np.maximum(1 - (arc(x, c0) / 0.9) ** 2, 0) ** 2
for nm, c0 in [("local_rare", np.pi), ("local_common", 0.0)]:
    val = integrate.quad(lambda x: (1 + 0.75 * np.cos(x)) / TWO_PI * 0.8 * bump(x, c0) / np.pi, 0, TWO_PI, points=[c0 % TWO_PI], limit=400)[0]
    rep["joint_tv_" + nm] = float(val)
rep["joint_tv_phase_eta08"] = 0.8 / np.pi
rep["attenuation_order4_order6"] = [float(np.exp(-k ** 2 * (0.35 ** 2 + 0.15 ** 2) / 2)) for k in (4, 6)]

json.dump(rep, open(os.path.join(out, "controle_analytique.json"), "w"), indent=2)
print(json.dumps(rep, indent=1))
