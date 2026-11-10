"""Tests of the crossing time and FCA functions on the exact model correlation, sampled on a lag
grid (no realizations):

    R(D, tau) = rho(tau) exp(-1/2 (D - v tau)^T W^-1 (D - v tau)),
    W = Rot(alpha) diag(l1^2, l2^2) Rot(alpha)^T,

plus the CorrelationDelegator on a two-pixel signal made of smoothed noise.
"""

import numpy as np
import pytest
from scipy.signal import fftconvolve

import velocity_estimation.full_correlation_analysis as fca
import velocity_estimation.utils as u

DT = 0.01
WINDOW = 10.0
LAGS = np.arange(-round(WINDOW / DT), round(WINDOW / DT) + 1) * DT

# A cubic spline through samples DT apart locates a crossing to about DT^4 times the fourth
# derivative of the difference of the correlations, 1e-8 here for velocities and widths of order
# one; the tolerance leaves two orders of magnitude.
TOL = 1e-5
# The shape fit is exact linear algebra on exact data.
SHAPE_TOL = 1e-9

# Separations of the pooled fit on a unit grid: the reference with its four neighbours, and the
# four pairs of a horizontal with a vertical neighbour.
NEIGHBOUR_SEPARATIONS = np.array([(1, 0), (-1, 0), (0, 1), (0, -1)], dtype=float)
DIAGONAL_SEPARATIONS = np.array([(-1, 1), (-1, -1), (1, 1), (1, -1)], dtype=float)
EIGHT_SEPARATIONS = np.vstack([NEIGHBOUR_SEPARATIONS, DIAGONAL_SEPARATIONS])
# One triangle: the reference with a right and an upper neighbour, and the pair of those two.
TRIANGLE_SEPARATIONS = np.array([(1, 0), (0, 1), (-1, 1)], dtype=float)


def rotation(alpha):
    return np.array([[np.cos(alpha), -np.sin(alpha)], [np.sin(alpha), np.cos(alpha)]])


def covariance(l1, l2, alpha):
    """W with the long axis l1 along (cos alpha, sin alpha)."""
    rot = rotation(alpha)
    return rot @ np.diag([l1**2, l2**2]) @ rot.T


def quadratic_form(separation, v, w_inv):
    """(D - v tau)^T W^-1 (D - v tau) on LAGS."""
    displacement = np.asarray(separation, dtype=float)[None, :] - np.outer(LAGS, v)
    return np.einsum("ij,jk,ik->i", displacement, w_inv, displacement)


def gaussian_correlation(separation, v, w, lifetime=np.inf):
    """R on LAGS for Gaussian pulses with the lifetime factor exp(-tau^2 / 2 lifetime^2)."""
    rho = np.exp(-(LAGS**2) / (2 * lifetime**2))
    return rho * np.exp(-quadratic_form(separation, v, np.linalg.inv(w)) / 2)


def non_gaussian_correlation(separation, v, width, lifetime):
    """R on LAGS for round pulses of shape (1 + |u|^2 / 4 width^2)^-2 with an exponential
    factor of the lag."""
    rho = np.exp(-np.abs(LAGS) / lifetime)
    return rho * (1 + quadratic_form(separation, v, np.eye(2) / width**2) / 4) ** -2


def zero_lag_correlations(separations, w):
    """R(D, 0) = exp(-1/2 D^T W^-1 D) for each separation."""
    w_inv = np.linalg.inv(w)
    return np.exp(-np.einsum("ij,jk,ik->i", separations, w_inv, separations) / 2)


def crossing_times(separations, v, w, lifetime=np.inf):
    r_auto = gaussian_correlation((0, 0), v, w, lifetime)
    return np.array(
        [
            fca.crossing_time(LAGS, r_auto, gaussian_correlation(d, v, w, lifetime))
            for d in separations
        ]
    )


def barber_pole_reading(e, v, w_inv):
    """What a crossing along the unit vector e reads, e^T W^-1 v / e^T W^-1 e."""
    return e @ w_inv @ v / (e @ w_inv @ e)


# 1. crossing_time ------------------------------------------------------------------------------


@pytest.mark.parametrize("v", [(1, 0.5), (-1, 0.5), (0.3, -1.2), (2, 2)])
@pytest.mark.parametrize("lifetime", [np.inf, 2.0, 0.5])
@pytest.mark.parametrize("separation", NEIGHBOUR_SEPARATIONS.tolist())
def test_crossing_time_round_gaussian_pulses(v, lifetime, separation):
    v, d = np.array(v, dtype=float), np.array(separation)
    w = np.eye(2)
    tau = fca.crossing_time(
        LAGS,
        gaussian_correlation((0, 0), v, w, lifetime),
        gaussian_correlation(d, v, w, lifetime),
    )
    component = d @ v / np.linalg.norm(d)
    assert np.isfinite(tau)
    assert np.sign(tau) == np.sign(component)
    assert np.linalg.norm(d) / (2 * tau) == pytest.approx(component, abs=TOL)


@pytest.mark.parametrize("v", [(1, 0.5), (-0.7, 0.4)])
@pytest.mark.parametrize("separation", NEIGHBOUR_SEPARATIONS.tolist())
def test_crossing_time_round_non_gaussian_pulses(v, separation):
    v, d = np.array(v, dtype=float), np.array(separation)
    width, lifetime = 1.0, 2.0
    tau = fca.crossing_time(
        LAGS,
        non_gaussian_correlation((0, 0), v, width, lifetime),
        non_gaussian_correlation(d, v, width, lifetime),
    )
    component = d @ v / np.linalg.norm(d)
    assert np.linalg.norm(d) / (2 * tau) == pytest.approx(component, abs=TOL)


def test_crossing_time_is_infinite_without_a_crossing_in_the_window():
    w, d = np.eye(2), np.array([1.0, 0.0])
    # Motion across the separation: R(D, tau) < R(0, tau) at every lag.
    v = np.array([0.0, 1.0])
    r_auto = gaussian_correlation((0, 0), v, w)
    assert fca.crossing_time(LAGS, r_auto, gaussian_correlation(d, v, w)) == np.inf
    # A small component along it: the crossing, at |D| / (2 e.v) = 50, lies beyond the window, so
    # the component reads zero with an error below |D| / (2 WINDOW).
    v = np.array([0.01, 1.0])
    r_auto = gaussian_correlation((0, 0), v, w)
    assert fca.crossing_time(LAGS, r_auto, gaussian_correlation(d, v, w)) == np.inf
    assert v[0] < np.linalg.norm(d) / (2 * WINDOW)


def test_crossing_time_is_nan_when_the_cross_correlation_is_not_below_the_autocorrelation():
    r_auto = np.exp(-(LAGS**2) / 2)
    assert np.isnan(fca.crossing_time(LAGS, r_auto, r_auto))
    assert np.isnan(fca.crossing_time(LAGS, r_auto, 1.1 * r_auto))
    assert np.isnan(fca.crossing_time(LAGS, r_auto, np.full_like(r_auto, np.nan)))


# 2. shape_from_zero_lag and shape_parameters ---------------------------------------------------


@pytest.mark.parametrize("aspect", [1, 1.5, 2, 4])
@pytest.mark.parametrize("tilt_deg", [0, 15, 30, 45, 60, 75, 90])
@pytest.mark.parametrize(
    "separations",
    [EIGHT_SEPARATIONS, TRIANGLE_SEPARATIONS],
    ids=["eight pairs", "one triangle"],
)
def test_shape_from_zero_lag(aspect, tilt_deg, separations):
    w = covariance(aspect, 1.0, np.radians(tilt_deg))
    w_inv = fca.shape_from_zero_lag(separations, zero_lag_correlations(separations, w))
    assert np.allclose(w_inv, np.linalg.inv(w), rtol=0, atol=SHAPE_TOL)
    l1, l2, tilt = fca.shape_parameters(w_inv)
    assert l1 == pytest.approx(aspect, abs=SHAPE_TOL)
    assert l2 == pytest.approx(1.0, abs=SHAPE_TOL)
    if aspect > 1:  # the tilt of a round shape is undefined
        assert tilt == pytest.approx(np.radians(tilt_deg), abs=1e-6)


@pytest.mark.parametrize("tilt_deg, expected_deg", [(120, -60), (-90, 90), (-30, -30)])
def test_shape_parameters_tilt_range(tilt_deg, expected_deg):
    w_inv = np.linalg.inv(covariance(2.0, 1.0, np.radians(tilt_deg)))
    assert fca.shape_parameters(w_inv)[2] == pytest.approx(
        np.radians(expected_deg), abs=1e-9
    )


# 3. velocity_from_crossings --------------------------------------------------------------------


@pytest.mark.parametrize(
    "aspect, tilt_deg",
    [(1, 0), (1.5, 45), (2, 30), (3, 75), (4, 0), (4, 30), (4, 60), (4, 90)],
)
@pytest.mark.parametrize("lifetime", [np.inf, 1.0, 0.3])
@pytest.mark.parametrize("v", [(-0.8, 0.5), (1, -0.6)])
def test_velocity_from_crossings(aspect, tilt_deg, lifetime, v):
    v = np.array(v, dtype=float)
    w = covariance(aspect, 1.0, np.radians(tilt_deg))
    w_inv_exact = np.linalg.inv(w)
    taus = crossing_times(NEIGHBOUR_SEPARATIONS, v, w, lifetime)
    assert np.all(
        np.isfinite(taus)
    )  # every crossing lies within the window in these cases

    # FCA: the crossings corrected by the measured shape give v for any lifetime.
    w_inv = fca.shape_from_zero_lag(
        EIGHT_SEPARATIONS, zero_lag_correlations(EIGHT_SEPARATIONS, w)
    )
    v_fca = fca.velocity_from_crossings(NEIGHBOUR_SEPARATIONS, taus, w_inv)
    assert np.allclose(v_fca, v, rtol=0, atol=TOL)

    # Crossing method: with W^-1 = I the fit returns the barber-pole reading along each axis.
    expected = [barber_pole_reading(e, v, w_inv_exact) for e in np.eye(2)]
    v_crossing = fca.velocity_from_crossings(NEIGHBOUR_SEPARATIONS, taus, np.eye(2))
    assert np.allclose(v_crossing, expected, rtol=0, atol=TOL)
    for d, tau in zip(NEIGHBOUR_SEPARATIONS, taus):
        e = d / np.linalg.norm(d)
        assert np.linalg.norm(d) / (2 * tau) == pytest.approx(
            barber_pole_reading(e, v, w_inv_exact), abs=TOL
        )


def test_velocity_from_crossings_infinite_crossing_reads_zero_component():
    v = np.array([0.0, 1.0])
    taus = crossing_times(NEIGHBOUR_SEPARATIONS, v, np.eye(2))
    assert np.isinf(taus[0]) and np.isinf(taus[1])
    v_crossing = fca.velocity_from_crossings(NEIGHBOUR_SEPARATIONS, taus, np.eye(2))
    assert np.allclose(v_crossing, v, rtol=0, atol=TOL)


# 4. spike_free_variance ------------------------------------------------------------------------


def test_spike_free_variance_removes_the_white_noise_spike():
    variance, noise_variance, correlation_time = 2.0, 0.7, 1.0
    acov = variance * np.exp(-(LAGS**2) / (2 * correlation_time**2))
    assert fca.spike_free_variance(acov, LAGS) == pytest.approx(variance, rel=1e-6)
    acov[np.argmin(np.abs(LAGS))] += noise_variance
    # The parabolic extrapolation from +-dt, +-2 dt misses the variance by dt^4 / 2 tau_c^4.
    assert fca.spike_free_variance(acov, LAGS) == pytest.approx(variance, rel=1e-6)


# 5. Invalid inputs -----------------------------------------------------------------------------


def test_shape_from_zero_lag_invalid_inputs():
    w = covariance(2.0, 1.0, 0.3)
    # Two directions only: no information on the off-diagonal element.
    assert (
        fca.shape_from_zero_lag(
            NEIGHBOUR_SEPARATIONS, zero_lag_correlations(NEIGHBOUR_SEPARATIONS, w)
        )
        is None
    )
    collinear = np.array([(1, 0), (2, 0), (3, 0)], dtype=float)
    assert (
        fca.shape_from_zero_lag(collinear, zero_lag_correlations(collinear, w)) is None
    )
    # Fewer than three pairs.
    two = TRIANGLE_SEPARATIONS[:2]
    assert fca.shape_from_zero_lag(two, zero_lag_correlations(two, w)) is None
    assert fca.shape_from_zero_lag([], []) is None
    # A correlation that is not positive.
    for bad in (0.0, -0.1, np.nan):
        r0 = zero_lag_correlations(EIGHT_SEPARATIONS, w)
        r0[0] = bad
        assert fca.shape_from_zero_lag(EIGHT_SEPARATIONS, r0) is None
    # Correlations in (0, 1) whose -2 ln R give Wxx = Wyy = 1, Wxy = 2: not positive definite.
    separations = np.array([(1, 0), (0, 1), (1, 1)], dtype=float)
    assert (
        fca.shape_from_zero_lag(separations, np.exp(-np.array([1, 1, 6]) / 2)) is None
    )


def test_velocity_from_crossings_invalid_inputs():
    along_x = np.array([(1, 0), (-1, 0)], dtype=float)
    assert fca.velocity_from_crossings(along_x, [0.5, -0.5], np.eye(2)) is None
    assert fca.velocity_from_crossings([], [], np.eye(2)) is None


# CorrelationDelegator --------------------------------------------------------------------------


class ArrayImagingData(u.ImagingDataInterface):
    """Signals (nx, ny, nt) on a unit grid."""

    def __init__(self, signals, dt):
        self.signals = signals
        self.dt = dt

    def get_shape(self):
        return self.signals.shape[:2]

    def get_signal(self, x, y):
        return self.signals[x, y]

    def get_dt(self):
        return self.dt

    def get_position(self, x, y):
        return float(x), float(y)

    def is_pixel_dead(self, x, y):
        return False


def test_correlation_delegator():
    dt, n, shift, width = 0.01, 100_000, 30, 20
    rng = np.random.default_rng(1)
    kernel = np.exp(-np.arange(-4 * width, 4 * width + 1) ** 2 / (2 * width**2))
    smooth = fftconvolve(rng.standard_normal(n + shift), kernel, mode="same")
    smooth /= smooth.std()
    reference, probe = (
        smooth[shift:],
        smooth[:-shift],
    )  # the probe sees it shift samples later
    signals = np.stack(
        [
            reference + 0.5 * rng.standard_normal(n),
            probe + 0.3 * rng.standard_normal(n),
        ]
    )[:, None, :]
    ds = ArrayImagingData(signals, dt)
    options = fca.FCAOptions(window=1.0)
    delegator = fca.CorrelationDelegator(options)

    lags, r_auto = delegator.autocorrelation((0, 0), ds)
    assert len(lags) == 201 and lags[0] == pytest.approx(-1.0) and lags[100] == 0
    assert r_auto[100] == 1
    # The spike-free variance is that of the noise-free signal, not of the noisy one (the noise's
    # own autocovariance at the first two lags fluctuates by 1 / sqrt(n) of its variance).
    assert delegator._pixel((0, 0), ds)[1] == pytest.approx(np.var(reference), rel=0.01)
    plain = fca.CorrelationDelegator(
        fca.FCAOptions(window=1.0, remove_noise_spike=False)
    )
    assert plain._pixel((0, 0), ds)[1] == pytest.approx(
        np.var(signals[0, 0]), rel=1e-12
    )
    assert plain.autocorrelation((0, 0), ds)[1][100] == 1

    # Structures move from pixel 0 to pixel 1, so R_10 peaks at a positive lag.
    lags, r_cross = delegator.cross_correlation((1, 0), (0, 0), ds)
    assert lags[np.argmax(r_cross)] == pytest.approx(shift * dt)
    assert np.max(r_cross) > 0.9
    assert delegator.zero_lag_correlation((1, 0), (0, 0), ds) == r_cross[100]

    # The reversed pair, served from the cache, equals a direct computation.
    _, r_reverse = delegator.cross_correlation((0, 0), (1, 0), ds)
    _, r_direct = fca.CorrelationDelegator(options, cache=False).cross_correlation(
        (0, 0), (1, 0), ds
    )
    assert np.allclose(r_reverse, r_direct)
    assert np.allclose(r_reverse, r_cross[::-1])


def test_zero_lag_without_full_correlation():
    """zero_lag_correlation computed directly, before any cross-correlation of the pair is
    cached, equals the zero lag of the full cross-correlation."""
    rng = np.random.default_rng(2)
    kernel = np.exp(-np.arange(-40, 41) ** 2 / (2 * 10**2))
    common = fftconvolve(rng.standard_normal(20_000), kernel, mode="same")
    signals = np.stack(
        [common + rng.standard_normal(20_000), common + rng.standard_normal(20_000)]
    )[:, None, :]
    ds = ArrayImagingData(signals, 0.01)
    direct = fca.CorrelationDelegator(fca.FCAOptions(window=0.5))
    r0 = direct.zero_lag_correlation((1, 0), (0, 0), ds)
    assert ((1, 0), (0, 0)) not in direct._pairs
    lags, r = fca.CorrelationDelegator(fca.FCAOptions(window=0.5)).cross_correlation(
        (1, 0), (0, 0), ds
    )
    assert r0 == pytest.approx(r[np.argmin(np.abs(lags))], rel=1e-10)
    # The reversed pair reads the cached value.
    assert direct.zero_lag_correlation((0, 0), (1, 0), ds) == r0


def test_non_positive_variance_gives_no_correlations():
    """A pixel whose autocovariance is negative at the first lag (the first difference of white
    noise: C(dt) = -C(0) / 2, C(2 dt) = 0) has a negative extrapolated variance. The delegator
    warns and leaves its correlations undefined, so the crossing time is NaN and the pixel
    gives no estimate instead of sign-flipped correlations."""
    rng = np.random.default_rng(3)
    white = rng.standard_normal(10_001)
    signals = np.stack([np.diff(white), rng.standard_normal(10_000)])[:, None, :]
    ds = ArrayImagingData(signals, 0.01)
    delegator = fca.CorrelationDelegator(fca.FCAOptions(window=0.2))
    with pytest.warns(UserWarning, match="non-positive variance"):
        lags, r_auto = delegator.autocorrelation((0, 0), ds)
    assert np.all(np.isnan(r_auto))
    _, r_cross = delegator.cross_correlation((1, 0), (0, 0), ds)
    assert np.all(np.isnan(r_cross))
    assert np.isnan(fca.crossing_time(lags, r_auto, r_cross))
    assert np.isnan(delegator.zero_lag_correlation((1, 0), (0, 0), ds))
