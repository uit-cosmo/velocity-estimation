"""Helpers shared by the realization tests: a blobmodel 2.x generator and the
closed-form readings of the estimators on its model."""

import numpy as np
from blobmodel import (
    BlobShapeEnum,
    BlobShapeImpl,
    DefaultBlobFactory,
    DistributionEnum,
    Geometry,
    Model,
)


def make_2d_realization(
    vx,
    vy,
    xpoints,
    ypoints,
    T=1000,
    K=1000,
    Lx=10,
    Ly=10,
    taup=np.inf,
    dt=0.01,
    t_lifetime=None,
    pad=None,
    wp=1,
    ws=1,
    theta=None,
    seed=None,
):
    """Realization of K Gaussian pulses of velocity (vx, vy) observed at the
    uniformly spaced pixels xpoints x ypoints, sampled at dt over [0, T).

    Returns the blobmodel dataset n(y, x, t) with the coordinates x and y equal
    to the requested points, as SyntheticBlobImagingDataInterface expects.

    The pulses have amplitude 1, widths wp along the direction theta from x
    (radians; None leaves them axis-aligned, which for wp == ws is round) and
    ws across it, the lifetime envelope exp(-((t - t_k) / t_lifetime)^2) when
    t_lifetime is given, and the drain time taup (inf: no draining). seed makes
    the realization reproducible.

    The view is the blobmodel Geometry grid. The process is homogeneous over it:

    - Along y the realization is periodic with period Ly, as with blobmodel 1.x:
      the pulse positions are seeded uniformly over a cell of length Ly and
      wrapped into it, with mirror copies at +-Ly, so every row sees the same
      pulse flux whatever vy. The cell is decoupled from the grid by setting
      the geometry's Ly after construction, since blobmodel ties the two and a
      cell as short as the view would make the mirror copies overlap the
      pulses. Ly must exceed the view's extent in y by 3 max(wp, ws) on each
      side, so that the periodic images change the correlations by less than
      exp(-9) ~ 1e-4; a smaller Ly raises. Without a lifetime the pulse density
      is K / (T Ly |vx|) per unit area, so the mean signal is K wp ws / (T Ly |vx|)
      and the variance is that over 2 pi (Campbell's theorem).
    - Along x the pulses are seeded (their position at t_k) uniformly over the
      view extended by pad on each side. With a finite lifetime a pulse is seen
      only within a few t_lifetime of t_k, so the seeding must reach a few
      |vx| t_lifetime, plus a width, beyond the view (blobmodel's
      docs/pulse_lifetime.rst): the default pad is
      4 (|vx| t_lifetime + max(wp, ws)). Without a lifetime every pulse exists
      at all times and any seeding gives the same homogeneous flux. Either way
      the first and last (Lx_view + 2 pad) / |vx| time units see fewer pulses,
      since none are seeded before t = 0 or after T: a transient of about 1 %
      of T = 1000 for the default view.

    Each pulse is summed only while it exceeds 1e-4 of its peak (blobmodel's
    speed_up truncation).

    Lx is unused: the x domain is the view. It stays in the signature so that
    positional calls written for the blobmodel 1.x helper keep working.
    """
    xpoints = np.asarray(xpoints, dtype=float)
    ypoints = np.asarray(ypoints, dtype=float)
    w_max = max(wp, ws)

    geometry = Geometry.from_arrays(
        xpoints, ypoints, np.arange(0, T, dt), periodic_y=True
    )
    if Ly < geometry.Ly + 6 * w_max:
        raise ValueError(
            f"Ly = {Ly} is too short a period for the view and the widths: it needs "
            f"at least {geometry.Ly + 6 * w_max:g}, the view's extent in y plus "
            f"3 max(wp, ws) on each side."
        )
    geometry.Ly = Ly

    if pad is None:
        pad = 4 * (abs(vx) * (t_lifetime or 0) + w_max)
    x_low, x_high = geometry.x0 - pad, geometry.x0 + geometry.Lx + pad

    bf = DefaultBlobFactory(t_drain=taup, t_lifetime=t_lifetime)
    bf.set_sampler("amplitude", DistributionEnum.deg, 1)
    bf.set_sampler("vx", DistributionEnum.deg, vx)
    bf.set_sampler("vy", DistributionEnum.deg, vy)
    bf.set_sampler("wp", DistributionEnum.deg, wp)
    bf.set_sampler("ws", DistributionEnum.deg, ws)
    bf.set_sampler("posx", lambda rng, n: rng.uniform(x_low, x_high, n))
    if theta is not None:
        bf.set_theta_setter(lambda: theta)

    model = Model(
        geometry=geometry,
        blob_shape=BlobShapeImpl(BlobShapeEnum.gaussian, BlobShapeEnum.gaussian),
        num_blobs=K,
        blob_factory=bf,
        verbose=False,
        seed=seed,
    )
    return model.make_realization(speed_up=True, truncation_error=1e-4)


def blob_w_inv(wp, ws, theta=0.0):
    """Inverse covariance W^-1 of the model cross-correlation
    R = rho exp(-1/2 (D - v tau)^T W^-1 (D - v tau)) for blobmodel's Gaussian
    pulses exp(-xb^2 / wp^2 - yb^2 / ws^2), xb along theta from x: their spatial
    autocorrelation is exp(-xb^2 / 2 wp^2 - yb^2 / 2 ws^2), so
    W = Rot(theta) diag(wp^2, ws^2) Rot(theta)^T."""
    rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    return rot @ np.diag([1 / wp**2, 1 / ws**2]) @ rot.T


def expected_tde_velocity(v, w_inv, t_lifetime=np.inf):
    """Reading of the three-point TDE on separations of equal length along x
    and y. With rho(tau) = exp(-tau^2 / 2 tau_p^2) the lag maximizing R(D, tau)
    is tau = v^T W^-1 D / (v^T W^-1 v + 1 / tau_p^2), and the 3TDE inverts the
    two delays as if they were the projections of v / |v|^2. Round pulses give
    v (1 + w^2 / (|v|^2 tau_p^2)); a frozen pattern (W^-1 v)(v^T W^-1 v) / |W^-1 v|^2.
    """
    v = np.asarray(v, dtype=float)
    w_inv_v = w_inv @ v
    return w_inv_v * (v @ w_inv_v + 1 / t_lifetime**2) / (w_inv_v @ w_inv_v)


def expected_crossing_velocity(v, w_inv):
    """Reading of the crossing times with W^-1 = I along x and y: the component
    along the unit vector e is e^T W^-1 v / e^T W^-1 e, for any lifetime. Round
    pulses give v."""
    v = np.asarray(v, dtype=float)
    return (w_inv @ v) / np.diag(w_inv)


def tilt_difference(a, b):
    """Difference of two axis angles in radians, modulo pi, in (-pi/2, pi/2]."""
    return (a - b + np.pi / 2) % np.pi - np.pi / 2
