"""The pulse model of the note (fca/note/main.tex) for the scans: a case of its parameters, its
cross-correlation function, and what each estimator gives on that function.

Units: the pulse size sqrt(l1 l2) = 1 and the mean speed |v| = 1, so times are in units of the
time a pulse takes to move its own size.
"""

from dataclasses import asdict, dataclass

import numpy as np

import velocity_estimation.full_correlation_analysis as fca

WINDOW = 8.0  # lag window of every estimate, [-WINDOW, WINDOW]
METHODS = ("2TDE", "3TDE", "CROSSING", "FCA")


@dataclass(frozen=True)
class Case:
    """One point of a scan.

    taup: pulse lifetime (inf: frozen pattern). aspect: l1 / l2 at l1 l2 = 1. tilt: angle of the
    long axis from x, degrees. angle: direction of the mean velocity from x, degrees. spread:
    standard deviation of each velocity component, normal and isotropic, over the mean speed.
    spacing: pixel separation Delta over the pulse size.
    """

    taup: float = np.inf
    aspect: float = 1.0
    tilt: float = 0.0
    angle: float = 0.0
    spread: float = 0.0
    spacing: float = 1.0

    def widths(self):
        return np.sqrt(self.aspect), 1 / np.sqrt(self.aspect)

    def velocity(self):
        beta = np.radians(self.angle)
        return np.array([np.cos(beta), np.sin(beta)])

    def as_dict(self):
        return asdict(self)


def rotation(degrees):
    c, s = np.cos(np.radians(degrees)), np.sin(np.radians(degrees))
    return np.array([[c, -s], [s, c]])


def correlation(case, separation, lags):
    """R(separation, lags) of the model: Gaussian pulses of covariance W, a Gaussian envelope of
    lifetime taup and normal velocities of mean v and covariance spread^2 I,

        R = rho(tau) sqrt(det W / det W(tau)) exp(-1/2 (D - v tau)^T W(tau)^-1 (D - v tau)),

    with rho = exp(-tau^2 / 2 taup^2) and W(tau) = W + spread^2 tau^2 I: the average of the
    single-velocity correlation over the velocities, which holds when the pulses are uniform in
    space and time whatever their velocity."""
    lags = np.asarray(lags, dtype=float)
    l1, l2 = case.widths()
    offset = np.asarray(separation, dtype=float)[:, None] - np.outer(case.velocity(), lags)
    along = rotation(case.tilt).T @ offset  # components along the long and the short axis
    variances = np.array([[l1**2], [l2**2]]) + (case.spread * lags) ** 2
    r = np.sqrt((l1 * l2) ** 2 / variances.prod(0)) * np.exp(-0.5 * (along**2 / variances).sum(0))
    if np.isfinite(case.taup):
        r = r * np.exp(-(lags**2) / (2 * case.taup**2))
    return r


def _delay(lags, r):
    """Lag of the maximum of r, refined on the parabola through its three highest samples."""
    k = int(np.argmax(r))
    if k in (0, len(r) - 1):
        return lags[k]
    a, b, c = r[k - 1 : k + 2]
    return lags[k] + 0.5 * (a - c) / (a - 2 * b + c) * (lags[1] - lags[0])


def expected(case, step=1e-3):
    """{method: (vx, vy)} that each estimator gives on the exact correlation function of the case,
    for the reference pixel with neighbours at the spacing along x and y: the definitions of the
    estimators, with the package's crossing and shape functions, applied to `correlation` on a
    fine lag grid. A two-point delay of zero leaves that 2TDE component undefined (NaN)."""
    n = int(round(WINDOW / step))
    lags = np.arange(-n, n + 1) * step
    d = case.spacing
    pairs = np.array([[d, 0.0], [0.0, d]])
    auto = correlation(case, (0, 0), lags)
    cross = [correlation(case, p, lags) for p in pairs]

    delays = np.array([_delay(lags, r) for r in cross])
    with np.errstate(divide="ignore"):
        two = np.where(np.abs(delays) > step / 10, d / delays, np.nan)
    three = d * delays / (delays**2).sum()

    crossings = [fca.crossing_time(lags, auto, r) for r in cross]
    crossing = fca.velocity_from_crossings(pairs, crossings, np.eye(2))

    shape_pairs = np.array([[d, 0], [0, d], [d, d], [d, -d]], dtype=float)
    r0 = [correlation(case, p, [0.0])[0] for p in shape_pairs]
    w_inv = fca.shape_from_zero_lag(shape_pairs, r0)
    full = fca.velocity_from_crossings(pairs, crossings, w_inv)

    return {"2TDE": two, "3TDE": three, "CROSSING": crossing, "FCA": full}
