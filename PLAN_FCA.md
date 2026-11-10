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

- **blobmodel 2.1.0** was released on 2026-11-10 (tag `v2.1.0` at `b9b4a52`, GitHub release, PyPI), so
  `blobmodel>=2.1.0` installs from PyPI.
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

## Results

Filled in as the steps complete (2026-11-10).

### Step 5: the realization tests

- blobmodel's Gaussian pulse is $\exp(-\theta^2)$, so its width parameter is the letter's $\ell$ directly ($\ell_1 = $
  `wp`, $\ell_2 = $ `ws`), and the lifetime of case 8 for $w = 1$, $\boldsymbol{v} = (1, 1)$ is
  $\tau_\mathrm{p} = 2\ell/|\boldsymbol{v}| = \sqrt{2}$. Measured 3TDE reading over 8 seeds: $1.252 \pm 0.012$ times
  $\boldsymbol{v}$ against the closed form 1.25; the crossing time and FCA read $\boldsymbol{v}$ within 0.05.
- Case 9 (aspect ratio 4, 30 degrees from $\boldsymbol{v} = (1, 1)$): keeping $R(\boldsymbol{\Delta}, 0)$ above
  `minimum_correlation` on the diagonal pairs needs $\ell_2 \geq 0.66$ pixels, hence $\ell_1 = 2.8$ pixels, and the
  FCA's $v_y$ is then ill-conditioned: the $y$ pair lies 15 degrees from the long axis ($R = 0.88$), its crossing is
  shallow, and $W$ amplifies the error of $\ln R$ by $\ell_1^2\sin^2\theta \approx 7$. Scatter over 8 seeds at
  $K = 2000$: 0.12 on $v_y$ (0.05 on $v_x$), largest error 0.27. The test bounds the FCA error by 0.4 (3.3 sigma; the
  barber-pole errors are 1.1 for the 3TDE and 1.9 for the crossing time) and checks $\ell_1/\ell_2$ within 20 % and the
  tilt within 0.1 rad. A second, well-conditioned tilted case (aspect ratio 2, $\ell_1 = 1$, $\ell_2 = 0.5$,
  $\theta = 30$ degrees, $\boldsymbol{v} = (1, 0)$) checks the 3TDE and the crossing time within 0.1 of their closed
  forms and the FCA within 0.1 of $\boldsymbol{v}$. The conditioning is now listed under the docs' limitations.
- Case 10: for noise variances up to the signal's, the velocity bias without the spike fix is at most 10 % (a neighbour
  with noise fraction $f$ against a noise-free reference moves the crossing to $R_\mathrm{cross}/R_\mathrm{auto} =
  (2 + f)/2\sqrt{1 + f}$), straddling the 0.1 tolerance (measured 0.08 to 0.12 over seeds), while $\ell_1$, $\ell_2$
  read 20 % low. The test asserts the ordering of the errors with and without the fix, the biased interval and the
  low widths; with the fix the FCA reads $\boldsymbol{v}$ within 0.04 and the widths within 5 %.
- Homogeneity along $y$: blobmodel ties the periodic cell to the grid, and a cell as short as the view would make the
  mirror copies overlap the pulses, so `make_2d_realization` sets the geometry's `Ly` (10 by default, at least the view
  plus three widths on each side) after construction, with the view as the grid. Measured per-pixel means and variances
  agree with the Poisson expectations within the statistical error (largest z-score 2.1 over 9 pixels and 10 seeds).
- `test_cond_av` passed in the baseline (blobmodel 1.0.1) and on 2.1; it is marginal (6 to 11 events) and fails in
  about 8 % of seeds. Not addressed.

### Step 6: the pooled fit against the mean of per-triangle solves

On the noisy realizations of case 10 (10 seeds, round pulses, $K = 1000$, $T = 1000$; the per-triangle FCA solves the
shape from the three pairs of one (horizontal, vertical) triangle and the velocity from its two crossings, exactly
determined, on the same correlations as the pooled fit):

| quantity | pooled: mean error, std over seeds | triangle mean: mean error, std | spread of the 4 triangles | pooled minus triangle mean |
|---|---|---|---|---|
| $v_x$ | 0.006, 0.020 | 0.006, 0.020 | 0.038 | $-0.0009 \pm 0.0015$ |
| $v_y$ | 0.008, 0.022 | 0.009, 0.023 | 0.039 | $-0.0014 \pm 0.0021$ |
| $\ell_1$ | 0.022, 0.018 | 0.033, 0.014 | 0.019 | $-0.011 \pm 0.009$ |
| $\ell_2$ | $-0.022$, 0.027 | $-0.027$, 0.022 | 0.017 | $0.004 \pm 0.008$ |

- The two agree to 0.001 per velocity component, 5 % of the scatter over seeds; neither is detectably biased in
  $\boldsymbol{v}$. With $W^{-1}$ near the identity the pooled least squares reduces to the crossing readings averaged
  over opposite neighbours, which is the triangle mean; they differ only at second order in $W^{-1}_{xy}$.
- The nonlinearity bias of averaging is real but small: the triangle mean's $\ell_1$ is 0.011 further above 1 and its
  $\ell_2$ 0.004 further below (the eigenvalues of four noisier $W^{-1}$ averaged after ordering); $\sqrt{\ell_1\ell_2}$
  is unbiased for both.
- Pooling does not change the noise bias: without the spike fix both read $v_x = 0.91$, $v_y = 0.92$ and
  $\ell_1, \ell_2 = 0.80, 0.77$; with it both are within 0.005 of their noise-free values seed by seed.
- No triangle failed in 120 (no NaN crossing, no pair below `minimum_correlation`, every $W^{-1}$ positive definite),
  so the pooled fit's robustness advantage is not exercised by this case.
- The pooled $W^{-1}_{xx}$, $W^{-1}_{yy}$ scatter 10 to 15 % more than the triangle mean's (not significant at 10 seeds)
  because the diagonal pairs, at $R = 0.37$, enter with the same weight as the axis pairs at $R = 0.61$; weighting the
  shape equations by $R$ would remove this. Not done.

Script and per-seed values: `scratchpad/step6/pooled_vs_triangles.py` of the session (not kept in the repository).

### Step 8: C-Mod shot 1160616027 through fusion_scripts

The branch (1.1.0, `64260fd`) is installed in fusion_scripts' venv in place of PyPI's 1.0.0 (restore with
`.venv/bin/pip install --no-deps --force-reinstall velocity_estimation==1.0.0`). The three fields were run with
`fields.load(1160616027)` and the hand-made dead-pixel mask (68 live pixels of 90; pitch 3.8 mm; $\mathrm{d}t$ = 0.5 µs,
0.3 s; autocorrelation e-folding lag 20 µs; the white-noise spike is 3 % of the variance at the median pixel, 36 % at
the worst). TDE as the twodca deck runs it (CC, 3-point, `minimum_cc_value` 0.5, running mean); CROSSING and FCA with
the default `FCAOptions`. Run times 78, 78 and 95 s; no warning from the crossing or shape fits.

| method | pixels | $|\boldsymbol{v}|$, m/s (median, IQR) | direction from radial (median, IQR) | median $v_r$, $v_z$ |
|---|---|---|---|---|
| TDE | 57 | 793 (739 to 992) | $-27.5^\circ$ ($-31.9$ to $-23.3$) | 721, $-389$ |
| CROSSING | 59 | 494 (437 to 550) | $-25.9^\circ$ ($-35.9$ to $-16.0$) | 442, $-232$ |
| FCA | 59 | 498 (365 to 581) | $+4.5^\circ$ ($-14.7$ to $+20.0$) | 486, $+28$ |

- At the 57 common pixels $|v_\mathrm{TDE}|/|v_\mathrm{CROSS}| = 1.76$ (1.39 to 2.10) at the same direction
  ($+1.3^\circ$): the decorrelation the letter predicts. The cross-correlation maxima lie at 1.5 to 3.5 µs, the crossing
  times at 5 to 30 µs, against the 20 µs autocorrelation time.
- The FCA shape: $\ell_1 = 0.89$ cm (2.4 pitches), $\ell_2 = 0.47$ cm (1.2 pitches), $\ell_1/\ell_2 = 1.9$, tilt
  $56^\circ$ from radial; the long axis is $88^\circ$ from the TDE and crossing vectors. The FCA turns the vector by
  $+28^\circ$ to nearly radial at the crossing's speed ($|v_\mathrm{CROSS}|/|v_\mathrm{FCA}| = 1.12$), flipping the
  median $v_z$ from $-232$ to $+28$ m/s. Its direction scatters more (IQR $35^\circ$ against $20^\circ$ for the
  crossing), as expected from $R(\boldsymbol{\Delta}, 0)$ of 0.76 to 0.88 at a pitch of 1.2 $\ell_2$ (the conditioning
  limitation above).
- Not in the plan, but decisive for reading the turn: the deck's cached 2DCA tracks on this shot, which assume neither a
  Gaussian shape nor a frozen pattern, read $v_r \approx 506$ to 510 and $v_z \approx +1$ to $+5$ m/s. The FCA agrees
  with them in direction ($+7^\circ$, IQR $-18$ to $+18$) and within about 15 % in speed; the TDE and crossing do not
  ($-27^\circ$, $-25^\circ$). Evidence for the FCA's correction, not proof: the 2DCA tracks large events only.
- Central pixel (x = 4, y = 5): $R(\boldsymbol{\Delta}, 0)$ 0.76 to 0.82 for the axis pairs, 0.55 to 0.88 for the
  diagonals; crossing times $-4.7$, $+5.4$, $+30$, $-15$ µs; $W^{-1}$ fit over 8 pairs with an RMS residual of 0.08 in
  $-2\ln R$; TDE (874, $-496$), CROSSING (390, $-109$), FCA (491, $+171$) m/s; the deck's 2DCA (707, $+145$).
- Follow-ups outside this step: `max_separation` = 2 for the shape fit's conditioning, and the other shots of
  `fields.SHOTS`. Fields, logs and the run script: `scratchpad/step8/` of the session.
- Installation note: fusion_scripts' `decorrelation/apd_check/simulation_column.py` imports
  `velocity_estimation.conditional_averaging`, which the repository dropped in `068bfd2` (2023) and only PyPI's 1.0.0
  wheel still shipped (from a stale `build/` directory); the clean 1.1.0 does not have it.

### The density scan: FCA fields of the nine shots

Through fusion_scripts' new `decorrelation/apd_check/fca_fields.py` (commit `3f97362` there; figures
`decorrelation/apd_check/figures/fca_field_<shot>.pdf` and `fca_fields_overview.pdf`, fields cached in
`decorrelation/cache/apd/apd<shot>_fca_field.nc`; regenerate with `python -m decorrelation.apd_check.fca_fields`).
Default `FCAOptions`, 54 to 59 of 68 live pixels estimated per shot, 45 to 128 s per shot.

| $f_\mathrm{GW}$ | shots | median $|\boldsymbol{v}|$, m/s | direction | $\ell_1$, $\ell_2$, cm | $\ell_1/\ell_2$ | tilt | confidence |
|---|---|---|---|---|---|---|---|
| 0.24 to 0.30 | 009, 011 | 250 to 290 | $-19^\circ$, $-63^\circ$ (scattered) | 0.6, 0.3 | 2.1 to 2.3 | $50^\circ$ | 0.63 |
| 0.44 to 0.46 | 016, 017, 018 | 290 to 300 | $+21^\circ$ to $+23^\circ$ | 0.65, 0.34 | 2.0 to 2.1 | $52^\circ$ | 0.65 |
| 0.65 | 022 | 430 | $+7^\circ$ | 0.85, 0.46 | 1.8 | $53^\circ$ | 0.79 |
| 0.80 | 025, 026, 027 | 500 to 515 | $+1^\circ$ to $+5^\circ$ | 0.9, 0.48 | 1.9 | $56^\circ$ | 0.82 |

- Speed, size and confidence rise with density; the aspect ratio (1.8 to 2.3) and the tilt of the long axis (50 to
  56 degrees from radial) do not change. The repeats at $f_\mathrm{GW}$ = 0.45 and 0.80 agree among themselves.
- The only failures are the shape fit at 1 to 5 edge pixels (one triangle) of the five low-density shots, where
  $\ell_2$ is below the pitch; none at $f_\mathrm{GW} \geq 0.65$. No crossing-time failure anywhere.
- Against the deck's cached fields: at $f_\mathrm{GW} \geq 0.65$ the step 8 pattern holds on every shot (FCA at 0.6 of
  the 3TDE speed, turned 28 to 37 degrees from it, within 15 % and 10 degrees of the 2DCA tracks). Below 0.65 the FCA
  reads 0.33 to 0.6 of the 3TDE speed and its direction agrees less with the 2DCA (30 degrees apart at 0.45), as
  expected where $\ell_2$ is below the pitch (the conditioning limitation).
- The two live columns nearest the separatrix carry the longest, mostly poloidal arrows, whose sign flips with
  density (downward at low density, upward at 0.80); they are also the low-confidence ones.

## Log

- 2026-11-10: plan written.
- 2026-11-10: the realization tests use blobmodel 2.1 (lifetime, start positions, tilt) instead of a generator of
  their own, at the author's request.
- 2026-11-10: blobmodel 2.1.0 released; installed from PyPI in a fresh venv, it builds the tilted, finite-lifetime case.
- 2026-11-10: steps 2 to 4 done (`a73d5ee`): the module, the delegator with the spike fix, the pipeline branch, tests 1 to
  6 (180 tests on the exact $R$).
- 2026-11-10: step 5 done (`269733f`): blobmodel 2.1.0 in the venv and as the `test` extra, `make_2d_realization` ported,
  tests 7 to 10 (the whole suite: 198 tests, 42 s). Deviations under Results.
- 2026-11-10: step 7 done (`e906294`): docs, README, version 1.1.0.
- 2026-11-10: review of the branch (eight findings; six fixed in `64260fd`: the variance guard, the direct zero-lag
  path, `FCAOptions` as a dataclass with `__all__`, the position lookup, the unused `Lx`; two declined: dropping
  zero-lag correlations above 1 would bias the shape fit, and falling back to the crossing where the shape fit fails
  would mix estimators in one field, so the hole is documented instead).
- 2026-11-10: step 6 done, results above: pooling and the triangle mean agree to 0.001 in $\boldsymbol{v}$ on case 10.
- 2026-11-10: step 8 done, results above: on shot 1160616027 the TDE reads 1.76 times the crossing speed, and the FCA
  turns the vector to the direction of the deck's 2DCA tracks. The branch is installed in fusion_scripts' venv.
- 2026-11-10: the `Sosnowsky/tripletail` branch (the LaTeX note `fca/note` and the manim slides `fca/presentation`)
  merged (`b4a9f49`), its worktree removed; the note gained a Limitations bullet on the conditioning of the shape fit
  with the spaced-antenna literature (`b612cc7`); the FCA fields of the nine density-scan shots drawn through
  fusion_scripts (results above).
