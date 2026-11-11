"""Realizations of the pulse model with blobmodel, estimated with the four methods of the package:
the scans of the note's section on realizations (fca/note/main.tex). From the repository root,

    venv/bin/python -m fca.scans.run [--workers 8] [--force] [scan ...]

runs the named scans (all by default) and writes fca/scans/results/<scan>.json, one row per
realization. A scan whose file exists is skipped unless --force is given. fca.scans.figures draws
the results.
"""

import argparse
import json
import platform
import subprocess
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from blobmodel import BlobShapeEnum, BlobShapeImpl, DefaultBlobFactory, DistributionEnum, Geometry, Model

import velocity_estimation.full_correlation_analysis as fca
import velocity_estimation.time_delay_estimation as tde
import velocity_estimation.two_dim_velocity_estimates as td
import velocity_estimation.utils as u
from fca.scans.model import WINDOW, Case

RESULTS = Path(__file__).parent / "results"

DT = 0.02
DURATION = 4000.0  # length of each record
LY = 30.0  # period along y: the view, WINDOW |v| and three widths on each side fit in it
FLUX = 0.1  # frozen pattern: pulses per unit time and unit length across the motion
DENSITY = 0.05  # finite lifetime: pulses per unit area within one lifetime, so that the pulses overlap as much in every case
SEEDS = 8  # records per case of a scan
ENSEMBLE = 240  # cases of the random ensemble, one record each

# The base of every scan but those that say otherwise: round frozen pulses moving along x with
# one velocity, pixels one pulse size apart.
ASPECTS_OF_TILT_SCAN = (2.0, 4.0)
SPREAD_LIFETIME = 5.0  # a spread needs a lifetime for the pulses to be uniform in space and time


def ensemble_cases(n=ENSEMBLE, seed=2026):
    """Cases with every parameter drawn at random: lifetime log-uniform in [1, 16], aspect ratio
    uniform in [1, 3], tilt uniform over the half turn, velocity direction uniform over the turn,
    spread uniform in [0, 0.3], spacing log-uniform in [0.5, 1.5]."""
    rng = np.random.default_rng(seed)
    return [
        Case(
            taup=float(np.exp(rng.uniform(np.log(1), np.log(16)))),
            aspect=float(rng.uniform(1, 3)),
            tilt=float(rng.uniform(-90, 90)),
            angle=float(rng.uniform(-180, 180)),
            spread=float(rng.uniform(0, 0.3)),
            spacing=float(np.exp(rng.uniform(np.log(0.5), np.log(1.5)))),
        )
        for _ in range(n)
    ]


SCANS = {
    # decorrelation: round pulses along x
    "lifetime": [Case(taup=t) for t in (0.5, 1, 2, 4, 8, 16, np.inf)],
    # misalignment: round frozen pulses
    "angle": [Case(angle=b) for b in (0, 15, 30, 45, 60, 75)],
    # barber pole: frozen elongated pulses along x
    "aspect": [Case(aspect=a, tilt=30) for a in (1, 1.5, 2, 3, 4)],
    "tilt": [Case(aspect=a, tilt=t) for a in ASPECTS_OF_TILT_SCAN for t in range(0, 91, 15)],
    # a spread of velocities, outside the model of the estimators
    "spread_round": [Case(taup=SPREAD_LIFETIME, spread=s) for s in (0, 0.1, 0.2, 0.3, 0.4, 0.5)],
    "spread_tilted": [
        Case(taup=SPREAD_LIFETIME, aspect=2, tilt=30, spread=s) for s in (0, 0.1, 0.2, 0.3, 0.4, 0.5)
    ],
    # pixel separation, with a lifetime and a tilt
    "spacing": [Case(taup=2, aspect=2, tilt=30, spacing=d) for d in (0.35, 0.5, 0.7, 1, 1.4, 2)],
    "ensemble": ensemble_cases(),
}


def realization(case, seed, duration=DURATION):
    """blobmodel dataset of the case on 3 x 3 pixels at the case's spacing.

    The pulses are Gaussian with unit amplitude. Along y the realization is periodic with period
    LY (blobmodel ties the period to the grid; it is set after construction, as in
    tests/test_utils.make_2d_realization). A frozen pattern has every pulse cross the view once,
    FLUX pulses per unit time and unit length. With a lifetime the pulses are seeded uniformly in
    time and over the view extended along x by three times the distance a fast pulse moves in a
    lifetime plus its width, DENSITY / taup pulses per unit area and time: uniform in space and time
    whatever their velocity, which the model correlation with a spread assumes."""
    l1, l2 = case.widths()
    vx, vy = case.velocity()
    points = case.spacing * np.arange(3.0)
    geometry = Geometry.from_arrays(points, points, np.arange(0, duration, DT), periodic_y=True)
    geometry.Ly = LY

    frozen = not np.isfinite(case.taup)
    if frozen:
        if case.spread > 0:
            raise ValueError("a velocity spread needs a finite lifetime")
        pad = 4 * l1
        k = FLUX * duration * LY
    else:
        pad = 3 * ((abs(vx) + 3 * case.spread) * case.taup + l1)
        k = DENSITY / case.taup * duration * LY * (geometry.Lx + 2 * pad)
    x_low, x_high = geometry.x0 - pad, geometry.x0 + geometry.Lx + pad

    bf = DefaultBlobFactory(t_drain=np.inf, t_lifetime=None if frozen else case.taup)
    bf.set_sampler("amplitude", DistributionEnum.deg, 1)
    bf.set_sampler("wp", DistributionEnum.deg, l1)
    bf.set_sampler("ws", DistributionEnum.deg, l2)
    bf.set_sampler("posx", lambda rng, n: rng.uniform(x_low, x_high, n))
    if case.spread > 0:
        bf.set_sampler("vx", lambda rng, n: rng.normal(vx, case.spread, n))
        bf.set_sampler("vy", lambda rng, n: rng.normal(vy, case.spread, n))
    else:
        bf.set_sampler("vx", DistributionEnum.deg, vx)
        bf.set_sampler("vy", DistributionEnum.deg, vy)
    tilt = np.radians(case.tilt)
    bf.set_theta_setter(lambda: tilt)

    model = Model(
        geometry=geometry,
        blob_shape=BlobShapeImpl(BlobShapeEnum.gaussian, BlobShapeEnum.gaussian),
        num_blobs=int(round(k)),
        blob_factory=bf,
        verbose=False,
        seed=seed,
    )
    return model.make_realization(speed_up=True, truncation_error=1e-4)


def estimate(ds):
    """{method: [vx, vy]} at the centre pixel, and the FCA shape [l1, l2, tilt in degrees]. The
    time delays are the maxima of the cross-correlations on a spline, with no lower bound on the
    maximum; every method uses the lag window [-WINDOW, WINDOW]. NaN where a method gives no
    estimate."""
    data = u.SyntheticBlobImagingDataInterface(ds)
    cc = tde.CCOptions(cc_window=WINDOW, minimum_cc_value=0, interpolate=True)
    options = {
        "2TDE": td.EstimationOptions(use_3point_method=False, cc_options=cc),
        "3TDE": td.EstimationOptions(use_3point_method=True, cc_options=cc),
        "CROSSING": td.EstimationOptions(
            velocity_method=td.VelocityMethod.CROSSING, fca_options=fca.FCAOptions(window=WINDOW)
        ),
        "FCA": td.EstimationOptions(
            velocity_method=td.VelocityMethod.FCA, fca_options=fca.FCAOptions(window=WINDOW)
        ),
    }
    delays = tde.TDEDelegator(tde.TDEMethod.CC, cc, True)
    correlations = fca.CorrelationDelegator(options["FCA"].fca_options, True)
    out, shape = {}, [np.nan] * 3
    for method, eo in options.items():
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pixel = td.estimate_velocities_for_pixel(1, 1, data, eo, delays, correlations)
        except ValueError:  # both delays zero
            pixel = td.PixelData()
        out[method] = [float(pixel.vx), float(pixel.vy)]
        if method == "FCA":
            shape = [float(pixel.l1), float(pixel.l2), float(np.degrees(pixel.tilt))]
    return out, shape


def run_one(task):
    case, seed = task
    start = time.time()
    estimates, shape = estimate(realization(case, seed))
    return {**case.as_dict(), "seed": seed, "estimates": estimates, "shape": shape, "seconds": round(time.time() - start, 1)}


def tasks(name):
    seeds = 1 if name == "ensemble" else SEEDS
    return [(case, 1000 * i + k) for i, case in enumerate(SCANS[name]) for k in range(seeds)]


def commit():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("scans", nargs="*", default=list(SCANS))
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    RESULTS.mkdir(exist_ok=True)
    with ProcessPoolExecutor(args.workers) as pool:
        for name in args.scans:
            path = RESULTS / f"{name}.json"
            if path.exists() and not args.force:
                print(f"{name}: {path} exists, skipped", flush=True)
                continue
            start = time.time()
            rows = list(pool.map(run_one, tasks(name)))
            settings = dict(dt=DT, duration=DURATION, Ly=LY, flux=FLUX, density=DENSITY, window=WINDOW)
            meta = dict(scan=name, commit=commit(), host=platform.node(), date=time.strftime("%Y-%m-%dT%H:%M:%S"), **settings)
            path.write_text(json.dumps({"meta": meta, "rows": rows}, indent=1))
            print(f"{name}: {len(rows)} realizations in {time.time() - start:.0f} s", flush=True)


if __name__ == "__main__":
    main()
