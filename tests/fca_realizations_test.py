"""Realization tests of the crossing time and FCA velocity estimates (PLAN_FCA.md, tests 7 to
10): blobmodel 2.1 realizations from tests/test_utils.make_2d_realization on a 3 x 3 view with
unit spacing, estimated at the reference pixel (1, 1) as in two_dim_velocity_estimates_test.py.
The expected readings of the biased estimators are the closed forms of tests/test_utils."""

import numpy as np

import test_utils as tu
import velocity_estimation.full_correlation_analysis as fca
import velocity_estimation.two_dim_velocity_estimates as td
import velocity_estimation.utils as u

X = Y = np.array([5.0, 6.0, 7.0])


def estimate(ds, method, window=None, remove_noise_spike=True, minimum_cc_value=None):
    """Estimate at the reference pixel (1, 1) with the given velocity method. window sets the
    lag window of both the time delays (cc_window) and the correlations of the crossings, and
    minimum_cc_value the smallest cross-correlation maximum the time delays accept; None keeps
    the defaults, 100 dt and 0.5."""
    options = td.EstimationOptions(
        velocity_method=method,
        fca_options=fca.FCAOptions(
            window=window, remove_noise_spike=remove_noise_spike
        ),
    )
    if window is not None:
        options.cc_options.cc_window = window
    if minimum_cc_value is not None:
        options.cc_options.minimum_cc_value = minimum_cc_value
    return td.estimate_velocities_for_pixel(
        1, 1, u.SyntheticBlobImagingDataInterface(ds), options
    )


def error(pixel_data, expected):
    """Largest error of the two velocity components."""
    return max(abs(pixel_data.vx - expected[0]), abs(pixel_data.vy - expected[1]))


def test_round_pulses():
    """Round pulses of width 1 moving with v = (1, 1), no lifetime: the crossing times and the
    FCA read v within the 0.1 tolerance of the 3TDE tests. Both scatter by 0.03 per component
    over seeds at K = 1000 (largest error over 8 seeds 0.044)."""
    v = (1, 1)
    ds = tu.make_2d_realization(*v, X, Y, seed=1)
    for method in (td.VelocityMethod.CROSSING, td.VelocityMethod.FCA):
        assert error(estimate(ds, method), v) < 0.1, method


def test_finite_lifetime():
    """Round pulses of width w = 1 with the lifetime tau_p = sqrt(2) = 2 w / |v|, v = (1, 1).
    The decorrelation pulls the maxima of the cross-correlations to shorter lags, so the 3TDE
    reads v (1 + w^2 / (|v|^2 tau_p^2)) = 1.25 v, outside the 0.1 tolerance of the frozen
    tests, while the crossing times and the FCA read v, the factor of the lag cancelling in the
    crossing. K = 8000 pulses over T = 4000: a pulse with a lifetime is seen only when seeded
    within about |v| tau_p of a pixel, so the pulses must be dense, and at this density the
    scatter is limited by the length of the series rather than by their number (it did not
    fall from K = 6000 to 16000 at T = 1000). Over 6 seeds the 3TDE scatters by 0.02 per
    component, the crossing by 0.03 and the FCA by 0.03 (largest errors 0.03, 0.05, 0.06).
    """
    v = (1, 1)
    t_lifetime = np.sqrt(2)
    ds = tu.make_2d_realization(*v, X, Y, T=4000, K=8000, t_lifetime=t_lifetime, seed=1)

    tde = estimate(ds, td.VelocityMethod.TDE)
    expected_tde = tu.expected_tde_velocity(v, tu.blob_w_inv(1, 1), t_lifetime)
    assert np.allclose(expected_tde, 1.25)
    assert error(tde, expected_tde) < 0.1
    assert error(tde, v) > 0.1

    assert error(estimate(ds, td.VelocityMethod.CROSSING), v) < 0.1
    assert error(estimate(ds, td.VelocityMethod.FCA), v) < 0.1


def test_tilted_elongated_pulses():
    """Frozen pattern of pulses of aspect ratio 4 (wp = 2.8, ws = 0.7) whose long axis is
    tilted 30 degrees from v = (1, 1), i.e. theta = 75 degrees from x. The barber pole turns
    the 3TDE towards the short axis and shortens it, (W^-1 v)(v^T W^-1 v) / |W^-1 v|^2 =
    (0.83, -0.13), and the crossing along each axis reads e^T W^-1 v / e^T W^-1 e = (0.75, -0.87);
    both match their closed forms (the 3TDE to 0.03 over 8 seeds at K = 2000). The FCA recovers
    the shape, l1 / l2 and the tilt, and removes the gross error of the velocity, but on this
    grid its vy is ill-conditioned: the y pair lies within 15 degrees of the long axis, so it
    is almost fully correlated (R(e_y, 0) = 0.88), its crossing is shallow (slope
    e_y^T W^-1 v = -0.22) and W amplifies the error by wp^2 sin^2(theta) ~ 7, a scatter of
    0.12 on vy (0.05 on vx) at K = 2000. The bound of 0.4 is 3.3 sigma (largest error over
    8 seeds 0.27), against the barber-pole errors of 1.1 (3TDE) and 1.9 (crossing). The x
    delay of the 3TDE, 1.18, needs a lag window wider than the default 100 dt = 1.0."""
    v = (1, 1)
    wp, ws, theta = 2.8, 0.7, np.deg2rad(75)
    ds = tu.make_2d_realization(
        *v, X, Y, K=2000, Ly=20, wp=wp, ws=ws, theta=theta, seed=1
    )
    w_inv = tu.blob_w_inv(wp, ws, theta)

    tde = estimate(ds, td.VelocityMethod.TDE, window=3.0)
    assert error(tde, tu.expected_tde_velocity(v, w_inv)) < 0.1
    assert error(tde, v) > 1.0

    crossing = estimate(ds, td.VelocityMethod.CROSSING, window=3.0)
    assert error(crossing, tu.expected_crossing_velocity(v, w_inv)) < 0.3
    assert error(crossing, v) > 1.0

    pixel = estimate(ds, td.VelocityMethod.FCA, window=3.0)
    assert error(pixel, v) < 0.4
    assert abs(pixel.l1 / pixel.l2 - wp / ws) < 0.2 * wp / ws
    assert abs(tu.tilt_difference(pixel.tilt, theta)) < 0.1


def test_tilted_pulses_well_conditioned():
    """Frozen pulses of aspect ratio 2 and of the size of the pixel spacing (wp = 1,
    ws = 0.5), long axis 30 degrees from v = (1, 0): a tilted case in which the FCA is well
    conditioned, both crossings having slopes e^T W^-1 v of order 1 (1.75 and -1.3) and W
    entries of order 1. The 3TDE reads (W^-1 v)(v^T W^-1 v) / |W^-1 v|^2 = (0.65, -0.48) and
    the crossing times e^T W^-1 v / e^T W^-1 e = (1, -0.40): the barber pole, a spurious
    velocity across the motion (the x crossing is exact for any shape when v is along x, the
    signals of two x-separated pixels being shifted copies). The FCA reads v within 0.1 and
    returns the aspect ratio within 20 % and the tilt within 0.1 rad. K = 2000: over 8 seeds
    the 3TDE scatters by 0.03 per component, the crossing by 0.03, the FCA by 0.02 (largest
    error 0.05), l1 / l2 by 0.09 and the tilt by 0.03 rad. The y crossing, at a lag of 1.25,
    and the 3TDE's x delay, 1.0, need a window wider than the default 100 dt = 1.0, and the
    cross-correlation across the short axis peaks at 0.32, below the 0.5 that the time delays
    require by default."""
    v = (1, 0)
    wp, ws, theta = 1.0, 0.5, np.deg2rad(30)
    ds = tu.make_2d_realization(*v, X, Y, K=2000, wp=wp, ws=ws, theta=theta, seed=1)
    w_inv = tu.blob_w_inv(wp, ws, theta)

    tde = estimate(ds, td.VelocityMethod.TDE, window=3.0, minimum_cc_value=0.2)
    assert error(tde, tu.expected_tde_velocity(v, w_inv)) < 0.1
    assert error(tde, v) > 0.3

    crossing = estimate(ds, td.VelocityMethod.CROSSING, window=3.0)
    assert error(crossing, tu.expected_crossing_velocity(v, w_inv)) < 0.1
    assert error(crossing, v) > 0.3

    pixel = estimate(ds, td.VelocityMethod.FCA, window=3.0)
    assert error(pixel, v) < 0.1
    assert abs(pixel.l1 / pixel.l2 - wp / ws) < 0.2 * wp / ws
    assert abs(tu.tilt_difference(pixel.tilt, theta)) < 0.1


def test_white_noise_spike():
    """Round pulses of width w = 1, v = (1, 1), plus white noise of a different variance at
    each pixel: none at the reference, 0.7 to 1.0 of the signal variance at its four
    neighbours and 0.3 to 0.6 at the corners. White noise adds a spike at zero lag to each
    pixel's autocovariance only. With the spike-free normalization (the default) the FCA reads
    v within 0.1 (largest error over 8 seeds 0.04) and l1, l2 within 10 % of w (5 %). Without
    it, a neighbour with noise fraction f against the noise-free reference moves their crossing
    to where R_cross / R_auto = (2 + f) / (2 sqrt(1 + f)), so the velocity reads low, 0.91 and
    0.92 for these fractions (the FCA measures 0.91 and 0.92 +- 0.015 over 8 seeds, the
    mean of its two components 0.89 to 0.94), and
    -2 ln R(D, 0) gains ln((1 + f_a)(1 + f_b)), so l1 and l2 read about 20 % low (0.80, 0.77).
    """
    v = (1, 1)
    # Noise variance over signal variance per pixel, indexed [y, x]; (1, 1) is the reference.
    fractions = np.array([[0.4, 0.9, 0.6], [1.0, 0.0, 0.8], [0.5, 0.7, 0.3]])
    ds = tu.make_2d_realization(*v, X, Y, seed=1)
    rng = np.random.default_rng(1)
    noise_std = np.sqrt(fractions * ds["n"].var("t").values)
    ds["n"] = ds["n"] + rng.normal(size=ds["n"].shape) * noise_std[:, :, np.newaxis]

    with_fix = estimate(ds, td.VelocityMethod.FCA)
    assert error(with_fix, v) < 0.1
    assert abs(with_fix.l1 - 1) < 0.1 and abs(with_fix.l2 - 1) < 0.1

    without_fix = estimate(ds, td.VelocityMethod.FCA, remove_noise_spike=False)
    assert error(without_fix, v) > error(with_fix, v)
    assert 0.85 < (without_fix.vx + without_fix.vy) / 2 < 0.96
    assert without_fix.l1 < 0.9 and without_fix.l2 < 0.9
