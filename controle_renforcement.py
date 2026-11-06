"""Independent CSV, kernel-CV and wind-score checks; does not call R."""
import csv
import numpy as np
from scipy.special import logsumexp, i0e, i1e
from scipy.optimize import brentq

def rows(name):
    with open(name, newline="") as stream:
        return list(csv.DictReader(stream))

for study in ["rotation", "sensibilite"]:
    raw = rows(f"renforcement/{study}_resultats.csv")
    summary = rows(f"renforcement/{study}_synthese.csv")
    grouping = ["condition", "model", "fidelity", "method"] if study == "rotation" else ["scenario", "n", "method"]
    groups = {}
    for row in raw:
        key = tuple(row[k] for k in grouping)
        groups.setdefault(key, []).append(row["reject"] == "TRUE")
    for row in summary:
        values = groups[tuple(row[k] for k in grouping)]
        assert len(values) == int(row["datasets"])
        assert sum(values) == int(row["rejections"])
        np.testing.assert_allclose(np.mean(values), float(row["rate"]), atol=1e-14, rtol=0)

pairs = rows("application/paires.csv")
days = sorted({row["date"] for row in pairs})
train = [row for row in pairs if row["date"] in days[:50]]
test = [row for row in pairs if row["date"] in days[50:]]
assert len(train) == 600 and len(test) == 276
x = np.array([float(row["x"]) for row in train])
y = np.array([float(row["y"]) for row in train])
day = np.array([row["date"] for row in train])
distance_x = x[:, None] - x
distance_y = y[:, None] - y

def log_kernel(distance, h):
    return logsumexp(np.stack([-0.5*((distance + 2*np.pi*j)/h)**2 - np.log(h*np.sqrt(2*np.pi)) for j in range(-2, 3)]), axis=0)

grid = rows("application/bandes_cv.csv")
values = []
for row in grid:
    lx = log_kernel(distance_x, float(row["hx"]))
    ly = log_kernel(distance_y, float(row["ht"]))
    lx[day[:, None] == day] = -np.inf
    value = (logsumexp(lx+ly, axis=1)-logsumexp(lx, axis=1)).sum()
    np.testing.assert_allclose(value, float(row["loglik"]), atol=2e-9, rtol=1e-12)
    values.append(value)
assert int(np.argmax(values)) == next(i for i, row in enumerate(grid) if row["selected"] == "TRUE")
mean_vector = np.exp(1j*(y-x)).mean()
kappa = brentq(lambda k: i1e(k)/i0e(k)-abs(mean_vector), 1e-8, 1000)
shift = np.angle(mean_vector)
fit = rows("application/persistance_ajustement.csv")[0]
np.testing.assert_allclose(kappa, float(fit["kappa"]), atol=1e-10)
score = rows("application/scores_paires.csv")
xt = np.array([float(row["x"]) for row in test])
yt = np.array([float(row["y"]) for row in test])
log_density = kappa*np.cos(yt-xt-shift)-np.log(2*np.pi)-np.log(i0e(kappa))-kappa
np.testing.assert_allclose(log_density, [float(row["log_persistence"]) for row in score], atol=1e-11)
for row in rows("renforcement/rotation_energie.csv"):
    rho, delta = float(row["rho"]), float(row["delta"])
    k = np.arange(1, 501)
    frequency = 2 if row["axial"] == "TRUE" else 1
    weight = rho**(2*k)*np.sin(frequency*k*delta)**2
    ratio = 1 if frequency == 2 else weight[k % 2 == 0].sum()/weight.sum()
    np.testing.assert_allclose(ratio, float(row["series"]), atol=1e-12)
print("Independent Python check passed: outcome counts, all 25 CV scores, persistence fit and scores, Fourier fractions.")
