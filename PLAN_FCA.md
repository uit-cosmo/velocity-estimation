# Plan: crossing time and full correlation analysis (FCA)

Branch `fca`, started 2026-11-10. Adds two velocity estimators to the per-pixel pipeline, next to the time-delay
estimation (TDE): the crossing time and FCA, as described in the letter by J. M. Losada and O. E. Garcia,
*Velocity estimation from three-point measurements: decorrelation and barber-pole biases and their correction by full
correlation analysis* (in preparation; working title).

## Why

For a superposition of pulses with a finite lifetime, elongated and tilted, the time delays are biased three ways:

- **Misalignment.** Two-point TDE overestimates the component along the separation when the motion is oblique to it.
  The 3TDE (the current default) removes this.
- **Decorrelation.** A finite lifetime, or any other factor of the lag alone, pulls the maximum of the
  cross-correlation towards shorter lags. Both the 2TDE and the 3TDE then read too high, and bringing the points
  closer does not help.
- **Barber pole.** Elongated structures tilted with respect to their motion turn the 3TDE towards their short axis
  and shorten it.

The crossing time removes the first two with two points; FCA removes all three with three points.

## The method

The model's cross-correlation, with $W$ the covariance of the pulse autocorrelation (widths $\ell_1, \ell_2$, tilt
$\alpha$):

$$R(\boldsymbol{\Delta}, \triangle_t) = \rho_\varphi(\triangle_t/\tau_\mathrm{p})\,
\exp\!\left[-\tfrac{1}{2}(\boldsymbol{\Delta} - \boldsymbol{v}\triangle_t)^\mathsf{T} W^{-1}
(\boldsymbol{\Delta} - \boldsymbol{v}\triangle_t)\right].$$

- **Crossing time.** $\tau_\mathrm{e}$ such that $R(0, \tau_\mathrm{e}) = R(\boldsymbol{\Delta}, \tau_\mathrm{e})$.
  The factor of the lag cancels, and
  $\boldsymbol{\Delta}^\mathsf{T} W^{-1} \boldsymbol{v} = \boldsymbol{\Delta}^\mathsf{T} W^{-1} \boldsymbol{\Delta} / 2\tau_\mathrm{e}$.
  For round pulses this is $\hat{\boldsymbol{e}}\cdot\boldsymbol{v} = \Delta/2\tau_\mathrm{e}$, exactly, for any lifetime
  and any pulse shape.
- **Shape.** At zero lag $-2\ln R(\boldsymbol{\Delta}, 0) = \boldsymbol{\Delta}^\mathsf{T} W^{-1} \boldsymbol{\Delta}$,
  linear in the three elements of $W^{-1}$.
- **FCA.** The crossing equations with the measured $W^{-1}$, solved for $\boldsymbol{v}$. With
  $W^{-1} = I$ they are the plain crossing time.

## Decisions (settled with the author, 2026-11-10)

1. **The crossing uses the mean of the two pixels' autocorrelations**, not the reference's alone. Swapping
   the two pixels then only changes the sign, so the cache can reuse reversed pairs, and an inhomogeneous field (the
   scrape-off layer's radial profile) does not favour one end of the pair. In the model the two choices agree.
2. **The crossing is the sampled one**: the first sign change of $R_\mathrm{auto} - R_\mathrm{cross}$ on each side of
   zero lag, refined on a cubic spline, the one nearest zero kept (port of `barberpole2d.crossing_lag` in
   fusion_scripts). Where there is none within the lag window, $1/\tau_\mathrm{e} = 0$: that component reads zero, an
   error of at most $\Delta/2T_\mathrm{window}$. A spurious crossing in the noisy tails is far from zero lag and so also
   reads a small component. The Gaussian fit near the peaks (G-INT, Holloway et al. 1997; Doviak et al. 2004), which
   the paper's barber-pole points use, assumes a Gaussian shape; it is added only if C-Mod data show the need.
3. **The shape is $-2\ln R(\boldsymbol{\Delta}, 0)$**, the letter's equation, which assumes Gaussian pulses. Briggs's
   calibration of the zero-lag values against the measured autocorrelation (exact for any elliptic frozen pattern) is a
   later extension.
4. **Spike-free normalization, on by default for the crossing time and FCA.** White noise adds a spike at zero lag to
   each pixel's autocovariance only. Each pixel's variance is replaced by its autocovariance extrapolated to zero lag,
   $C^*(0) = [4C(\mathrm{d}t) - C(2\,\mathrm{d}t)]/3$ (port of `moments2d.condition_r`), and every correlation is
   normalized with it: $R_{ab}(\tau) = C_{ab}(\tau)/\sqrt{C^*_{aa}(0)\,C^*_{bb}(0)}$. Without it, uncorrelated noise
   lowers every $R(\boldsymbol{\Delta}, 0)$ and adds a constant to every $-2\ln R$, which biases $W^{-1}$, and unequal
   noise between two pixels shifts their crossing. TDE is unchanged.
5. **The measured shape is returned**: $\ell_1$, $\ell_2$ (major, minor) and the tilt $\alpha$ of the long axis from the
   first coordinate, from the eigenvectors of $W^{-1}$. New `PixelData` fields, NaN unless the method is FCA.
6. **One pooled fit per reference pixel**, not one solve per triangle averaged (the 3TDE's way). See below.

## Pooling over the neighbours

The current pipeline takes the horizontal and the vertical neighbours of the reference pixel (`_find_neighbors`, with
`NeighbourOptions`), solves the 3TDE for each (horizontal, vertical) triangle and averages the velocities. FCA instead
collects the same pairs and solves once:

- **Shape**, least squares for $(W^{-1}_{xx}, W^{-1}_{xy}, W^{-1}_{yy})$: one equation
  $\Delta_x^2 W^{-1}_{xx} + 2\Delta_x\Delta_y W^{-1}_{xy} + \Delta_y^2 W^{-1}_{yy} = -2\ln R(\boldsymbol{\Delta}, 0)$
  per pair of the reference with each neighbour, and of each horizontal neighbour with each vertical one. Up to eight
  equations for three unknowns. Pairs of two horizontal (or two vertical) neighbours, at twice the separation, are left
  out: for non-Gaussian pulses $-2\ln R$ there measures a different curvature.
- **Velocity**, least squares for $\boldsymbol{v}$: one crossing equation per pair of the reference with a neighbour.
  Up to four equations for two unknowns.

Why: with both diagonals, $W^{-1}_{xy} = [\ln R(\Delta, -\Delta) - \ln R(\Delta, \Delta)]/2\Delta^2$ comes from two
correlations, while one triangle gets it from three logarithms; triangles that fail the validity checks on their own
would otherwise be dropped selectively; and the mean of four nonlinear FCA solutions is not the FCA of the data. At a
boundary pixel with one triangle the same code is exactly determined. Opposite neighbours keep the estimate centred on
the reference pixel. Arbitrary separations (the tilted C-Mod grid, `max_separation > 1`) need nothing special.

Both directions are required: with only horizontal or only vertical neighbours there is no estimate, as for the 3TDE.

## API

Backward compatible: `EstimationOptions()` keeps its default (3TDE), and `eo.use_3point_method = ...` set after
construction, as imaging_methods and fusion_scripts do, keeps working.

```python
class VelocityMethod(Enum):
    TDE = 1       # time delays; 2- or 3-point by use_3point_method (current behaviour)
    CROSSING = 2  # crossing times, W^-1 = I: the component along each separation
    FCA = 3       # crossing times corrected by the measured shape

class FCAOptions:
    window: float = None             # lags in [-window, window]; None: 100 dt, as CCOptions.cc_window
    remove_noise_spike: bool = True  # decision 4
    minimum_correlation: float = 0.1 # pairs with R(D, 0) below it are left out of the shape fit
```

- `EstimationOptions` gets `velocity_method` (default `VelocityMethod.TDE`) and `fca_options`, both used by
  `CROSSING` and `FCA`; `method` and the TDE options are ignored there. Defaults are built inside `__init__`, not as
  shared default instances (imaging_methods works around that sharing for `CCOptions` today).
- `PixelData` gets `l1`, `l2`, `tilt` (NaN by default); `MovieData` gets `get_l1`, `get_l2`, `get_tilt`.
- `confidence` is the smallest maximum of the reference's cross-correlations, as for the CC TDE; `events` stays 0.
- Sign convention as in the TDE: $R(\boldsymbol{\Delta}, \tau)$ correlates the pixel at
  $\boldsymbol{\Delta} = \boldsymbol{r}_1 - \boldsymbol{r}_0$ at $t + \tau$ with the reference at $t$, i.e.
  `corr_fun(signal(p1), signal(p0))`, which peaks at $\tau > 0$ when structures move from $p_0$ to $p_1$.

## Code

New module `velocity_estimation/full_correlation_analysis.py`. Pure functions of correlation arrays, so that they can
be tested on an exact $R$:

- `spike_free_variance(acov, lags)`: $C^*(0)$ from the autocovariance at $\pm\mathrm{d}t, \pm 2\,\mathrm{d}t$.
- `crossing_time(lags, r_auto, r_cross)`: decision 2; `np.inf` when there is no crossing.
- `shape_from_zero_lag(separations, r0)`: least squares for $W^{-1}$; `None` if fewer than three non-collinear
  separations, any $R \le 0$ among those used, or $W^{-1}$ not positive definite.
- `velocity_from_crossings(separations, crossing_times, w_inv)`: least squares for $\boldsymbol{v}$; `None` if the
  separations do not span the plane.
- `shape_parameters(w_inv)`: $(\ell_1, \ell_2, \alpha)$.

A `CorrelationDelegator`, like `TDEDelegator`: caches each pixel's spike-free variance and autocorrelation and each
pair's cross-covariance over the window, reusing a reversed pair as $R_{ba}(\tau) = R_{ab}(-\tau)$.

In `two_dim_velocity_estimates.py`, `estimate_velocities_for_pixel` branches on `velocity_method`: neighbour selection
is shared; `CROSSING` and `FCA` build the pairs of the pooled fit and call the functions above. Invalid inputs give no
estimate and a warning, as the TDE does.

## Tests

Exact $R$ (fast, no realizations):

1. `crossing_time` on the Gaussian $R$ sampled at a finite `dt`: $\Delta/2\tau_\mathrm{e} = \hat{\boldsymbol{e}}\cdot\boldsymbol{v}$
   for round pulses, oblique motion, finite and infinite lifetime, both signs; also exact for a round non-Gaussian
   shape, e.g. $(1 + |\boldsymbol{u}|^2/4)^{-2}$ with an exponential factor of the lag.
2. `shape_from_zero_lag` recovers $W^{-1}$ for aspect ratios 1 to 4 and tilts 0 to 90 degrees, from eight pairs and from
   one triangle.
3. `velocity_from_crossings` with the measured $W^{-1}$ gives $\boldsymbol{v}$ for the same shapes and any lifetime;
   with $W^{-1} = I$ it gives the crossings' barber-pole reading,
   $\hat{\boldsymbol{e}}^\mathsf{T} W^{-1}\boldsymbol{v}/\hat{\boldsymbol{e}}^\mathsf{T} W^{-1}\hat{\boldsymbol{e}}$ per separation.
4. `spike_free_variance` on an autocovariance plus a white-noise spike returns the noise-free variance.
5. Invalid inputs (collinear separations, $R(\boldsymbol{\Delta}, 0) \le 0$, $W^{-1}$ not positive definite) give `None`.
6. Backward compatibility: the default options, `use_3point_method` set after construction, the existing tests.

Realizations (the pipeline end to end) with blobmodel 2.1, which has what the letter's model needs:

- `DefaultBlobFactory(t_lifetime=...)`: the Gaussian envelope $\exp[-((t - t_\mathrm{init})/\tau_\mathrm{p})^2]$, the
  letter's $\varphi$ with $\tau_\mathrm{p}$ = `t_lifetime`;
- `set_sampler("posx", lambda rng, n: rng.uniform(x0 - pad, x0 + Lx + pad, n))`: start positions across a domain wider
  than the view by a few times $v\tau_\mathrm{p}$, so that the process is homogeneous over the view (blobmodel's own
  `docs/pulse_lifetime.rst`);
- `set_theta_setter(lambda: alpha)` with widths `wp`, `ws`: elliptic pulses tilted from $x$ independently of
  $\boldsymbol{v}$;
- `periodic_y=True`, so that the view is homogeneous along $y$ too (to check when porting).

The repository's blobmodel 1.0.1 cannot produce these cases: it starts every pulse at $x = 0$ and decays it with the
time since then, so all pulses reach a given pixel scaled by the same factor and the normalized $R$ shows no
decorrelation; and it tilts pulses only along their velocity. Checked 2026-11-10 with the `main` checkout: a 5 x 5 view,
4000 pulses tilted 30 degrees with $\tau_\mathrm{p} = 2$, 40 000 samples, in 2 s.

`tests/test_utils.make_2d_realization` is ported to the 2.x API (`Geometry`, `set_sampler`,
`BlobShapeImpl(BlobShapeEnum.gaussian, BlobShapeEnum.gaussian)`), with parameters for the lifetime, the start-position
pad, the widths and the tilt. Its `update_geometry`, which patches 1.x internals to place the pixels, goes: the view is a
`Geometry` grid (`x0`, `y0`, `Nx`, `Ny`). The existing tests keep their cases and tolerances. Cases:

7. Round pulses, $\boldsymbol{v} = (1, 1)$: `CROSSING` and `FCA` within the tolerance of the current 3TDE tests.
8. A finite lifetime, $\tau_\mathrm{p} = 2\ell/v$ with $\ell$ the blob width: the 3TDE reads high by about
   $1 + \ell^2/v^2\tau_\mathrm{p}^2$;
   `CROSSING` and `FCA` within tolerance.
9. A frozen pattern of elongated pulses tilted 30 degrees from $\boldsymbol{v}$, aspect ratio 4: the 3TDE and
   `CROSSING` off by their closed forms, `FCA` within tolerance, and the returned $\ell_1/\ell_2$ and tilt close to the
   inputs.
10. White noise of different variance at each pixel (up to the signal's): `FCA` within tolerance with the spike fix, and
    measurably biased without it, so that the test shows the fix is needed.

One-off check, results recorded below: the pooled fit against the mean of per-triangle FCA solves, scatter over seeds,
on the noisy realizations of case 10.

## Dependencies

- **blobmodel 2.1.0 is not released.** It is on `main` of uit-cosmo/2d-propagating-blobs (pyproject version 2.1.0),
  but the last tag, GitHub release and PyPI version are 2.0.0. Release it before `fca` is merged; until then the
  repository's venv installs the checkout (`pip install -e ~/Git/2d-propagating-blobs`).
- blobmodel is imported only by the tests, not by the package. It moves from the runtime `dependencies` in
  `pyproject.toml` to a `test` extra (`blobmodel>=2.1.0`, `pytest`), so that imaging_methods and fusion_scripts do not
  need 2.1.0 to install velocity_estimation. `requirements.txt` gets `blobmodel>=2.1.0`, together with the author's
  uncommitted change to it.

## Docs

- `docs/estimation_options.rst`: `VelocityMethod` and `FCAOptions` tables.
- `docs/usage.rst`: choosing the method; reading `get_l1`, `get_l2`, `get_tilt`.
- `docs/index.rst` and `README.md`: the three biases and the two estimators, citing the letter, Briggs (1950), Phillips
  and Spencer (1955) and Briggs (1984).
- Limitations: $-2\ln R$ assumes Gaussian pulses; the spike fix removes white noise only, not coherent pickup (C-Mod's
  lines at 200 to 1000 kHz) or cross-talk between neighbouring channels, which raises $R(\boldsymbol{\Delta}, 0)$ and the
  cross-correlations; with a velocity spread the shape changes with the lag, $W(\tau_\mathrm{e}) = W + \Sigma_V
  \tau_\mathrm{e}^2$, which the zero-lag shape does not contain; the extrapolated zero lag needs a sampling time well
  below the correlation time.
- Version 1.1.0.

## Steps

1. Branch and plan (this file).
2. `full_correlation_analysis.py` with tests 1 to 5.
3. `CorrelationDelegator` with the spike fix.
4. Options, `PixelData` fields and the pipeline branch; test 6.
5. blobmodel 2.1 in the venv and the `test` extra; `make_2d_realization` ported, the existing tests passing on it;
   tests 7 to 10.
6. The pooled-against-triangles check; results below.
7. Docs and version.
8. Outside the repository: the TDE, crossing and FCA fields on a C-Mod shot through fusion_scripts, whose venv
   installs velocity_estimation from site-packages, so the branch has to be installed there first.

## Out of scope

- Pooling the 3TDE (its delays are linear in the slowness $\boldsymbol{v}/|\boldsymbol{v}|^2$, so it could be pooled
  the same way, but that changes existing results).
- G-INT, Briggs's autocorrelation calibration of the shape, and the paper's lag-dependent shape with a velocity spread,
  $R(\boldsymbol{\Delta}, \tau)R(\boldsymbol{\Delta}, -\tau)/R(0, \tau)^2$.
- A crossing time from conditional averages.

## State at the start

- On `main` (`b3b6945`) 12 tests pass and `test_cond_av` already fails; not addressed here.
- `requirements.txt` has an uncommitted change of the author's (pins relaxed to minimums), carried onto this branch
  and left uncommitted.

## Log

- 2026-11-10: plan written.
- 2026-11-10: the realization tests use blobmodel 2.1 (lifetime, start positions, tilt) instead of a generator of
  their own, at the author's request.
