"""Crossing times and full correlation analysis (FCA).

For a superposition of pulses moving with velocity v, with W the covariance of the pulse
autocorrelation (widths l1 >= l2, tilt alpha of the long axis from the first coordinate) and a
lifetime factor rho of the lag alone, the cross-correlation between the reference pixel and a pixel
at separation D is

    R(D, tau) = rho(tau / tau_p) exp(-1/2 (D - v tau)^T W^-1 (D - v tau)).

The crossing time tau_e of a pair, R(0, tau_e) = R(D, tau_e), does not depend on rho and gives
D^T W^-1 v = D^T W^-1 D / (2 tau_e); the zero-lag correlations give W^-1 through
-2 ln R(D, 0) = D^T W^-1 D. The functions below take correlation arrays, so that they can be tested
on the exact R; CorrelationDelegator computes and caches the correlations of a dataset.

Sign convention as in the time delay estimation: R(D, tau) with D = r1 - r0 correlates the pixel
at r1 at time t + tau with the reference at r0 at time t, so it peaks at tau > 0 when structures
move from r0 to r1.
"""

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import brentq

from velocity_estimation import correlation as cf
from velocity_estimation import utils


@dataclass
class FCAOptions:
    def __init__(
        self,
        window: float = None,
        remove_noise_spike: bool = True,
        minimum_correlation: float = 0.1,
    ):
        """Options of the crossing time and FCA velocity estimates.

        - window: float, the correlations are computed for time lags in [-window, window]. If None,
        100 sampling times are used, as in CCOptions.
        - remove_noise_spike: bool, if True each pixel's variance is replaced by its autocovariance
        extrapolated to zero lag, [4 C(dt) - C(2 dt)] / 3, which removes the spike that white noise
        adds at zero lag, and every correlation is normalized with it.
        - minimum_correlation: float, pairs whose zero-lag correlation is below it are left out of
        the shape fit of the FCA.
        """
        self.window = window
        self.remove_noise_spike = remove_noise_spike
        self.minimum_correlation = minimum_correlation

    def __str__(self):
        """Return a string representation of the FCAOptions object."""
        return (
            f"Window: {self.window}, Remove Noise Spike: {self.remove_noise_spike}, "
            f"Minimum Correlation: {self.minimum_correlation}"
        )


def _zero_lag_index(lags):
    return int(np.argmin(np.abs(lags)))


def spike_free_variance(acov, lags):
    """Autocovariance extrapolated to zero lag, [4 C(dt) - C(2 dt)] / 3 with C(k dt) the mean of
    acov at +-k dt: the variance without the spike that white noise adds at zero lag. lags is
    uniformly spaced and contains 0."""
    mid = _zero_lag_index(lags)
    c1 = (acov[mid - 1] + acov[mid + 1]) / 2
    c2 = (acov[mid - 2] + acov[mid + 2]) / 2
    return float((4 * c1 - c2) / 3)


def crossing_time(lags, r_auto, r_cross):
    """Lag nearest zero, of either sign, at which the cross-correlation r_cross reaches the
    autocorrelation r_auto: the first sign change of r_cross - r_auto on each side of zero lag,
    refined on a cubic spline through the samples. np.inf when there is no crossing within the
    lags (the velocity component along the separation then reads zero), np.nan when r_cross is
    not below r_auto at zero lag."""
    diff = np.asarray(r_cross, dtype=float) - np.asarray(r_auto, dtype=float)
    mid = _zero_lag_index(lags)
    if not diff[mid] < 0:
        return np.nan
    spline = InterpolatedUnivariateSpline(lags, diff, k=3)
    roots = []
    for step in (1, -1):
        k = mid
        while 0 <= k + step < len(lags) and diff[k + step] < 0:
            k += step
        if not 0 <= k + step < len(lags):
            continue
        if diff[k + step] == 0:
            roots.append(lags[k + step])
        else:
            roots.append(brentq(spline, *sorted((lags[k], lags[k + step]))))
    return float(min(roots, key=abs)) if roots else np.inf


def shape_from_zero_lag(separations, r0):
    """Least-squares W^-1 from the zero-lag correlations r0 of pairs at the given separations
    (n, 2): one equation Dx^2 Wxx + 2 Dx Dy Wxy + Dy^2 Wyy = -2 ln r0 per pair. Returns the
    symmetric 2 x 2 matrix, or None with fewer than three pairs, a correlation that is not
    positive, fewer than three distinct directions among the separations, or a W^-1 that is not
    positive definite."""
    separations = np.asarray(separations, dtype=float).reshape(-1, 2)
    r0 = np.asarray(r0, dtype=float)
    if len(r0) < 3 or not np.all(r0 > 0):
        return None
    dx, dy = separations[:, 0], separations[:, 1]
    design = np.column_stack([dx**2, 2 * dx * dy, dy**2])
    if np.linalg.matrix_rank(design) < 3:
        return None
    wxx, wxy, wyy = np.linalg.lstsq(design, -2 * np.log(r0), rcond=None)[0]
    w_inv = np.array([[wxx, wxy], [wxy, wyy]])
    if np.any(np.linalg.eigvalsh(w_inv) <= 0):
        return None
    return w_inv


def velocity_from_crossings(separations, crossing_times, w_inv):
    """Least-squares velocity from the crossing times of pairs at the given separations (n, 2):
    one equation D^T W^-1 v = D^T W^-1 D / (2 tau) per pair, an infinite crossing time giving a
    zero right-hand side. With w_inv the identity this is the crossing time method, which reads
    the component of the velocity along each separation. Returns v (2,), or None if the
    separations do not span the plane."""
    separations = np.asarray(separations, dtype=float).reshape(-1, 2)
    crossing_times = np.asarray(crossing_times, dtype=float)
    rows = separations @ w_inv
    if len(rows) < 2 or np.linalg.matrix_rank(rows) < 2:
        return None
    rhs = np.sum(rows * separations, axis=1) / (2 * crossing_times)
    return np.linalg.lstsq(rows, rhs, rcond=None)[0]


def shape_parameters(w_inv):
    """Widths l1 >= l2 and tilt of the long axis from the first coordinate, in radians in
    (-pi/2, pi/2], from the eigenvectors of W^-1: l = 1 / sqrt(eigenvalue), the long axis along
    the eigenvector of the smallest."""
    eigenvalues, eigenvectors = np.linalg.eigh(w_inv)
    l1, l2 = 1 / np.sqrt(eigenvalues)
    ex, ey = eigenvectors[:, 0]
    tilt = np.pi / 2 - (np.pi / 2 - np.arctan2(ey, ex)) % np.pi  # an axis: modulo pi
    return float(l1), float(l2), float(tilt)


class CorrelationDelegator:
    """Computes the correlations of a dataset for the crossing time and FCA estimates over the lag
    window of the options, normalized with each pixel's spike-free variance when the options ask
    for it. When cache is True, keeps each pixel's variance and autocovariance and each pair's
    cross-covariance, reusing a reversed pair as R_ba(tau) = R_ab(-tau)."""

    def __init__(self, options: FCAOptions = None, cache: bool = True):
        self.options = options if options is not None else FCAOptions()
        self.cache = cache
        self._pixels = {}
        self._pairs = {}

    def _window(self, lags, dt):
        window = self.options.window if self.options.window is not None else 100 * dt
        return np.abs(lags) <= window

    def _pixel(self, p, ds: utils.ImagingDataInterface):
        """(lags, C*(0), autocovariance over the window with C*(0) at zero lag) of pixel p."""
        key = tuple(p)
        if key in self._pixels:
            return self._pixels[key]
        dt = ds.get_dt()
        signal = ds.get_signal(p[0], p[1])
        signal = signal - signal.mean()
        lags, acov = cf.corr_fun(signal, signal, dt, norm=False, biased=True)
        if self.options.remove_noise_spike:
            variance = spike_free_variance(acov, lags)
        else:
            variance = float(acov[_zero_lag_index(lags)])
        window = self._window(lags, dt)
        lags, acov = lags[window], acov[window]
        acov[_zero_lag_index(lags)] = variance
        result = (lags, variance, acov)
        if self.cache:
            self._pixels[key] = result
        return result

    def _cross_covariance(self, p1, p0, ds: utils.ImagingDataInterface):
        """Covariance of pixel p1 at t + tau with pixel p0 at t over the window."""
        key, reverse = (tuple(p1), tuple(p0)), (tuple(p0), tuple(p1))
        if key in self._pairs:
            return self._pairs[key]
        if reverse in self._pairs:
            return self._pairs[reverse][::-1]
        dt = ds.get_dt()
        x = ds.get_signal(p1[0], p1[1])
        y = ds.get_signal(p0[0], p0[1])
        lags, ccov = cf.corr_fun(
            x - x.mean(), y - y.mean(), dt, norm=False, biased=True
        )
        ccov = ccov[self._window(lags, dt)]
        if self.cache:
            self._pairs[key] = ccov
        return ccov

    def autocorrelation(self, p, ds: utils.ImagingDataInterface):
        """(lags, R) of pixel p with itself over the window, with R(0) = 1."""
        lags, variance, acov = self._pixel(p, ds)
        return lags, acov / variance

    def cross_correlation(self, p1, p0, ds: utils.ImagingDataInterface):
        """(lags, R) of pixel p1 at t + tau with pixel p0 at t over the window, normalized with
        the two pixels' variances; it peaks at tau > 0 when structures move from p0 to p1.
        """
        lags, variance1, _ = self._pixel(p1, ds)
        _, variance0, _ = self._pixel(p0, ds)
        ccov = self._cross_covariance(p1, p0, ds)
        return lags, ccov / np.sqrt(variance1 * variance0)

    def zero_lag_correlation(self, p1, p0, ds: utils.ImagingDataInterface):
        """R(D, 0) of the pair."""
        lags, r = self.cross_correlation(p1, p0, ds)
        return float(r[_zero_lag_index(lags)])
