# Crossing time and full correlation analysis: handout

State as of 2026-11-10, branch `fca` at `b4b653c`: twelve commits on `main` (`b3b6945`), not pushed. This file says what
the branch adds to velocity_estimation, why it took this shape, what it gave on synthetic and C-Mod data, where the
code, tests, documents and figures are, and what is left. The decisions and the full result tables are in
`PLAN_FCA.md`; the derivations are in the note `fca/note/main.tex`; the package docs in `docs/` are the user-facing
reference.

## What the branch adds

Version 1.1.0 adds two velocity estimators to the per-pixel pipeline of `two_dim_velocity_estimates.py`, next to the
time-delay estimation (TDE), which stays the default. Both work on the cross-correlations $R(\boldsymbol\Delta, \tau)$
of the reference pixel with its horizontal and vertical neighbours (the 3TDE's neighbours, from `NeighbourOptions`),
for the model of the letter by Losada and Garcia (in preparation): a superposition of pulses with velocity
$\boldsymbol v$, lifetime $\tau_\mathrm p$ and a Gaussian spatial autocorrelation of covariance $W$ (widths
$\ell_1 \geq \ell_2$, tilt $\alpha$ of the long axis from the first coordinate),

$$R(\boldsymbol\Delta, \tau) = \rho(\tau/\tau_\mathrm p)\,\exp\!\left[-\tfrac12(\boldsymbol\Delta - \boldsymbol v\tau)^\mathsf T W^{-1}(\boldsymbol\Delta - \boldsymbol v\tau)\right].$$

| method (`VelocityMethod`) | uses | gives | biased by |
|---|---|---|---|
| `TDE`, two points | the lag of the maximum of $R$ of one pair | $\Delta/\tau_\mathrm{max}$ along the separation | misalignment, decorrelation, barber pole |
| `TDE`, three points (default) | the two lags of a (horizontal, vertical) triangle, averaged over the triangles | $\boldsymbol v$ | decorrelation, barber pole |
| `CROSSING` | the crossing time of each pair, $R(0, \tau_\mathrm e) = R(\boldsymbol\Delta, \tau_\mathrm e)$, which gives $\boldsymbol\Delta\cdot\boldsymbol v = \Delta^2/2\tau_\mathrm e$ | $\boldsymbol v$, least squares over the pairs | barber pole |
| `FCA` | the crossing times and the zero-lag correlations, $-2\ln R(\boldsymbol\Delta, 0) = \boldsymbol\Delta^\mathsf T W^{-1}\boldsymbol\Delta$ | $\boldsymbol v$, $\ell_1$, $\ell_2$, $\alpha$ | nothing in the model; the conditioning, under Limitations in the docs and the note |

The biases. *Misalignment*: a two-point delay overestimates the component along the separation when the motion is
oblique to it. *Decorrelation*: a finite lifetime, or any factor of the lag alone, pulls the maximum of $R$ to shorter
lags, so every TDE gives too high a speed, $\boldsymbol v(1 + w^2/|\boldsymbol v|^2\tau_\mathrm p^2)$ for round Gaussian
pulses of width $w$, and bringing the points closer does not help. *Barber pole*: elongated pulses tilted with respect
to their motion turn the 3TDE towards their short axis, to
$(W^{-1}\boldsymbol v)(\boldsymbol v^\mathsf T W^{-1}\boldsymbol v)/|W^{-1}\boldsymbol v|^2$ for a frozen pattern. At the
crossing time the lag factor $\rho$ cancels, which removes the decorrelation bias exactly, for any pulse shape and
lifetime; the measured $W^{-1}$ removes the barber pole.

Also on the branch: the LaTeX note `fca/note/main.tex` (revtex4-2: the model, the crossing time, the shape from the
zero-lag correlations, FCA, the pooled fit on the grid, the correlations from measured signals, usage, limitations) and
manim slides of it in `fca/presentation/`, both from the author's `Sosnowsky/tripletail` branch, merged in `b4a9f49`.
The note gained a Limitations bullet on the conditioning of the shape fit, with the spaced-antenna literature
(`b612cc7`; see Literature below).

## Why this shape (decisions with the author, 2026-11-10)

- **The crossing uses the mean of the two pixels' autocorrelations.** Swapping the pixels then only changes the sign,
  so the cache reuses reversed pairs, and a radial profile does not favour one end of the pair.
- **The crossing is the sampled one**: the first sign change of $R_\mathrm{cross} - R_\mathrm{auto}$ on each side of zero
  lag, refined on a cubic spline, the nearest to zero kept (port of fusion_scripts' `barberpole2d.crossing_lag`). No
  crossing within the window gives $1/\tau_\mathrm e = 0$: that component counts as zero, an error of at most
  $\Delta/2T_\mathrm{window}$. No Gaussian fit near the peaks (G-INT) unless C-Mod data show the need.
- **The shape is $-2\ln R(\boldsymbol\Delta, 0)$**, which assumes Gaussian pulses. Briggs's calibration of the zero-lag
  values against the autocorrelation is a later extension.
- **Spike-free normalization, on by default.** White noise adds a spike at zero lag to each pixel's autocovariance
  only. Each pixel's variance is replaced by $C^*(0) = [4C(\mathrm dt) - C(2\,\mathrm dt)]/3$ (port of
  `moments2d.condition_r`) and every correlation is normalized with it. Without it, noise adds a constant to every
  $-2\ln R$, which biases $W^{-1}$, and unequal noise between two pixels shifts their crossing. The TDE is unchanged.
- **One pooled least-squares fit per reference pixel**, not one solve per triangle averaged as the 3TDE does: up to
  eight shape equations (the reference with each neighbour, each horizontal neighbour with each vertical one) for the
  three elements of $W^{-1}$, and up to four crossing equations for $\boldsymbol v$. Pairs of two horizontal (or two
  vertical) neighbours are left out: at twice the separation a non-Gaussian $-2\ln R$ measures a different curvature.
  Checked against the triangle mean (Results).
- **The measured shape is returned**: `PixelData.l1`, `l2`, `tilt` (radians from the first coordinate, in
  $(-\pi/2, \pi/2]$), NaN unless the method is FCA.
- **Backward compatible.** `EstimationOptions()` keeps the 3TDE; `eo.use_3point_method = ...` set after construction
  keeps working (imaging_methods and fusion_scripts do that); option defaults are built per instance, not shared.
- **blobmodel 2.1 for the realization tests**, at the author's request, instead of a generator of our own: it has the
  lifetime envelope, start positions over a wider domain and a tilt independent of the velocity, which 1.0.1 lacked
  (it started every pulse at $x = 0$ and decayed it since launch, which does not decorrelate). blobmodel moves from
  the runtime dependencies to the `test` extra, so imaging_methods and fusion_scripts install 1.1.0 without it.
- **Declined in review**: dropping zero-lag correlations above 1 (noise lets the spike-free normalization exceed 1,
  and dropping only the upward fluctuations would bias the fit), and falling back to the crossing time where the shape
  fit fails (it would mix two estimators in one field; the hole is documented instead).

## Results

### Tests

201 pass in about 55 s (`venv/bin/python -m pytest tests -q -p no:cacheprovider`). `test_cond_av` is marginal (6 to
11 events) and fails in about 8 % of seeds, as before the branch; not addressed.

- `tests/full_correlation_analysis_test.py`: the functions on the exact $R$ sampled at a finite $\mathrm dt$: crossing
  times for round Gaussian and non-Gaussian pulses, both signs, finite and infinite lifetime; $W^{-1}$ for aspect
  ratios 1 to 4 and tilts 0 to 90 degrees, from eight pairs and from one triangle; the FCA velocity and the
  crossing's barber-pole closed form; the spike removal; invalid inputs; the delegator's sign convention, cache and
  variance guard.
- `tests/estimation_options_test.py`: the defaults, options set after construction, fresh instances, `PixelData` NaNs,
  the package exports.
- `tests/fca_realizations_test.py`, blobmodel 2.1 realizations through the whole pipeline:

| case | input | checks |
|---|---|---|
| round pulses | $\boldsymbol v = (1, 1)$, $w = 1$ | 3TDE, crossing and FCA within 0.1 of $\boldsymbol v$ |
| finite lifetime | $\tau_\mathrm p = \sqrt 2 = 2\ell/\lvert\boldsymbol v\rvert$, $T = 4000$, $K = 8000$ | 3TDE within 0.1 of $1.25\,\boldsymbol v$ and more than 0.1 from $\boldsymbol v$; crossing and FCA within 0.1 of $\boldsymbol v$ |
| tilted, aspect ratio 4 | $\ell_1 = 2.8$, $\ell_2 = 0.7$, $75^\circ$ from $x$ ($30^\circ$ from $\boldsymbol v$), frozen | 3TDE and crossing off by their closed forms (errors 1.1 and 1.9); FCA within 0.4; $\ell_1/\ell_2$ within 20 %, tilt within 0.1 rad |
| tilted, aspect ratio 2 | $\ell_1 = 1$, $\ell_2 = 0.5$, $30^\circ$, $\boldsymbol v = (1, 0)$ | 3TDE and crossing within 0.1 of their closed forms; FCA within 0.1 of $\boldsymbol v$ |
| white noise | variance 0 at the reference, 0.7 to 1 of the signal's at the neighbours | with the spike fix FCA within 0.1 and widths within 0.1; without it a larger error, components in 0.85 to 0.96, widths below 0.9 |

Deviations from the plan, recorded in its Results: the aspect-ratio-4 case cannot meet 0.1 (ill-conditioned: the $y$
pair lies $15^\circ$ from the long axis at $R = 0.88$, and $W$ amplifies the error of $\ln R$ by about 7; scatter 0.12
over seeds), hence the 0.4 bound (3.3 sigma) and the added well-conditioned case; and the noise case's velocity bias
without the fix straddles 0.1 (closed form $(2 + f)/2\sqrt{1 + f}$ for a neighbour with noise fraction $f$), so the
test asserts the ordering of the errors and the biased interval instead of a threshold.

### The pooled fit against the mean of per-triangle solves

On the noisy realizations, 10 seeds: the two agree to 0.001 per velocity component, 5 % of the scatter over seeds, and
neither is biased in $\boldsymbol v$; the triangle mean's $\ell_1$ is 0.011 higher and its $\ell_2$ 0.004 lower
(ordered eigenvalues averaged); pooling does not change the noise bias (both give $v \approx 0.91$ without the spike
fix); no triangle failed, so the pooled fit's robustness is not exercised there. The pooled $W^{-1}_{xx}$,
$W^{-1}_{yy}$ scatter 10 to 15 % more, because the diagonal pairs ($R = 0.37$) weigh as much as the axis pairs
($R = 0.61$); weighting the shape equations by $R$ would remove this and is not done.

### C-Mod shot 1160616027

Through fusion_scripts: 68 live pixels of 90, pixel spacing 3.8 mm, $\mathrm dt$ = 0.5 µs over 0.3 s, autocorrelation
time 20 µs; the white-noise spike is 3 % of the variance at the median pixel and 36 % at the worst. TDE as the deck
runs it (CC, three points, `minimum_cc_value` 0.5), crossing and FCA at the default `FCAOptions`; about 80 s per field.

| method | pixels | $\lvert\boldsymbol v\rvert$ in m/s, median (IQR) | direction from radial, median (IQR) |
|---|---|---|---|
| 3TDE | 57 | 793 (739 to 992) | $-27.5^\circ$ ($-31.9$ to $-23.3$) |
| CROSSING | 59 | 494 (437 to 550) | $-25.9^\circ$ ($-35.9$ to $-16.0$) |
| FCA | 59 | 498 (365 to 581) | $+4.5^\circ$ ($-14.7$ to $+20.0$) |

- The 3TDE gives 1.76 times the crossing speed at the same direction: the decorrelation the letter predicts (the maxima
  lie at 1.5 to 3.5 µs, the crossings at 5 to 30 µs).
- The FCA shape is $\ell_1 = 0.89$ cm, $\ell_2 = 0.47$ cm (1.2 pixel spacings), tilt $56^\circ$ from radial, the long
  axis $88^\circ$ from the TDE vector. The FCA turns the vector by $28^\circ$ to nearly radial at the crossing's speed,
  and agrees in direction, and within about 15 % in speed, with the deck's cached 2DCA tracks, which assume neither a
  Gaussian shape nor a frozen pattern; the TDE and the crossing do not. Evidence for the correction, not proof: the 2DCA
  tracks large events only.
- The FCA direction scatters more (IQR $35^\circ$ against $20^\circ$), as expected with $R(\boldsymbol\Delta, 0)$ of
  0.76 to 0.88 at this spacing.

### The density scan

fusion_scripts' `fca_fields.py`, nine shots at the default options, 54 to 59 of 68 live pixels estimated per shot,
45 to 128 s per shot:

| $f_\mathrm{GW}$ | shots | median $\lvert\boldsymbol v\rvert$, m/s | direction | $\ell_1$, $\ell_2$, cm | $\ell_1/\ell_2$ | tilt | confidence |
|---|---|---|---|---|---|---|---|
| 0.24 to 0.30 | 009, 011 | 250 to 290 | $-19^\circ$, $-63^\circ$ (scattered) | 0.6, 0.3 | 2.1 to 2.3 | $50^\circ$ | 0.63 |
| 0.44 to 0.46 | 016, 017, 018 | 290 to 300 | $+21^\circ$ to $+23^\circ$ | 0.65, 0.34 | 2.0 to 2.1 | $52^\circ$ | 0.65 |
| 0.65 | 022 | 430 | $+7^\circ$ | 0.85, 0.46 | 1.8 | $53^\circ$ | 0.79 |
| 0.80 | 025, 026, 027 | 500 to 515 | $+1^\circ$ to $+5^\circ$ | 0.9, 0.48 | 1.9 | $56^\circ$ | 0.82 |

- Speed, size and confidence rise with density; the aspect ratio and the tilt of the long axis do not. The repeats at
  the same density agree among themselves.
- At $f_\mathrm{GW} \geq 0.65$ the pattern of shot 027 holds on every shot (FCA at 0.6 of the 3TDE speed, turned 28 to
  37 degrees from it, within 15 % and 10 degrees of the 2DCA tracks). Below it $\ell_2$ is below the pixel spacing, the
  FCA gives 0.33 to 0.6 of the 3TDE speed and its direction agrees less with the 2DCA ($30^\circ$ apart at 0.45): the
  conditioning limitation.
- The only failures are the shape fit at 1 to 5 edge pixels of the low-density shots; no crossing-time failure
  anywhere. The two live columns nearest the separatrix carry the longest, mostly poloidal, low-confidence arrows,
  whose sign flips with density.

## Where things are

| Path | What |
|---|---|
| `velocity_estimation/full_correlation_analysis.py` | `FCAOptions`; `spike_free_variance`, `crossing_time`, `shape_from_zero_lag`, `velocity_from_crossings`, `shape_parameters`, pure functions of correlation arrays; `CorrelationDelegator`, the cached correlations of a dataset. |
| `velocity_estimation/two_dim_velocity_estimates.py` | `VelocityMethod`; `EstimationOptions.velocity_method` and `.fca_options`; `PixelData.l1`, `.l2`, `.tilt`; `MovieData.get_l1`, `get_l2`, `get_tilt`; `_estimate_velocities_from_crossings` (the pooled fit) and the branch in `estimate_velocities_for_pixel`. |
| `tests/test_utils.py` | `make_2d_realization` on blobmodel 2.1 (lifetime, pad, widths, tilt, seed; the periodic cell `Ly` decoupled from the view) and the closed forms `blob_w_inv`, `expected_tde_velocity`, `expected_crossing_velocity`, `tilt_difference`. |
| `tests/full_correlation_analysis_test.py`, `estimation_options_test.py`, `fca_realizations_test.py` | The new tests. `two_dim_velocity_estimates_test.py` keeps the old ones on the ported generator. |
| `docs/index.rst`, `usage.rst`, `estimation_options.rst`, `README.md` | The biases and the two estimators, choosing the method, the shape getters and the width convention, the option tables, the limitations, the references. |
| `PLAN_FCA.md` | The plan: decisions, the full result tables (steps 5, 6, 8 and the scan) and the log. |
| `fca/note/main.tex`, `bibliography.bib` | The note. `fca/presentation/slides.py` and `README.md`: the manim slides and how to render them. |
| fusion_scripts `decorrelation/apd_check/fca_fields.py` | `compute(shot)`, cached in `decorrelation/cache/apd/apd<shot>_fca_field.nc` (vr, vz, confidence, l1, l2, tilt, dead); `summary`, `table`, `fig_field`, `fig_overview`, `main`. Commit `3f97362` on `main`. |
| fusion_scripts `decorrelation/apd_check/figures/fca_field_<shot>.pdf`, `fca_fields_overview.pdf` | The nine fields and the 3 x 3 overview, with their `provenance.json` entries. |
| fusion_scripts `decorrelation/apd_check/fields.py`, `twodca_manuscript/tde_field.py` | The preprocessed APD records (`fields.load(shot)`, `fields.DEAD`, `fields.SHOTS`), the deck's TDE settings and its cached 2DCA tracks. |

The step 6 and step 8 scripts and the literature record lived in the session's scratchpad and are gone with it; their
results are in the plan and below.

## How to run

```bash
# velocity-estimation, from the repo root
venv/bin/pip install -e '.[test]'                        # blobmodel>=2.1.0 and pytest
venv/bin/python -m pytest tests -q -p no:cacheprovider   # 201 tests, about 55 s
venv/bin/python -m sphinx -b html docs <out>             # the docs; the autosummary warning about lumache predates the branch
cd fca/note && latexmk -pdf main                         # the note
cd fca/presentation && .venv/bin/manim-slides render --quality=h slides.py <scenes>   # the slides, see its README

# fusion_scripts, from its root; its venv installs velocity_estimation from site-packages, not editable
rm -rf /home/sosno/Git/velocity-estimation/build         # or setuptools packs stale files into the wheel
.venv/bin/pip install --no-deps --force-reinstall /home/sosno/Git/velocity-estimation
.venv/bin/python -m decorrelation.apd_check.fca_fields              # every shot: computes what is not cached, draws all
.venv/bin/python -m decorrelation.apd_check.fca_fields 1160616027   # one shot
.venv/bin/pip install --no-deps --force-reinstall velocity_estimation==1.0.0   # back to the PyPI release
```

Usage, with `import velocity_estimation as ve`:

```python
eo = ve.EstimationOptions(velocity_method=ve.VelocityMethod.FCA)  # or CROSSING; the TDE options are then ignored
eo.fca_options.window = None              # lags in [-window, window]; None: 100 dt, as CCOptions.cc_window
eo.fca_options.remove_noise_spike = True
eo.fca_options.minimum_correlation = 0.1  # pairs with R(D, 0) below it stay out of the shape fit
movie_data = ve.estimate_velocity_field(ds, eo)
vx, vy = movie_data.get_vx(), movie_data.get_vy()
l1, l2, tilt = movie_data.get_l1(), movie_data.get_l2(), movie_data.get_tilt()
```

## Conventions and gotchas

- **Widths.** For a correlation $\exp(-x^2/2\sigma^2)$, $\ell = \sqrt2\,\sigma$. blobmodel's Gaussian pulse is
  $\exp(-\theta^2)$, so its `wp`, `ws` are $\ell_1$, $\ell_2$ directly, and `t_lifetime` is $\tau_\mathrm p$, with
  $\rho(\tau) = \exp(-\tau^2/2\tau_\mathrm p^2)$. `t_drain` does not decorrelate.
- **Sign.** $R(\boldsymbol\Delta, \tau)$ with $\boldsymbol\Delta = \boldsymbol r_1 - \boldsymbol r_0$ is
  `corr_fun(signal(p1), signal(p0))`, which peaks at $\tau > 0$ when structures move from $p_0$ to $p_1$, as in the TDE.
- **Window.** The default 100 dt is the TDE's. On realizations of elongated tilted pulses the maxima lie outside it and
  the 3TDE gives nothing, so those tests pass `window=3.0` to both `cc_options` and `fca_options`. `minimum_correlation`
  0.1 keeps the pairs across the short axis out of the shape fit when $\ell_2$ is below the pixel spacing.
- **Holes.** A pixel without a horizontal or a vertical neighbour, a NaN crossing (the cross-correlation not below the
  autocorrelation at zero lag, or an undefined variance) or a failed shape fit (fewer than three directions, $R \le 0$,
  $W^{-1}$ not positive definite) gives no estimate and a warning, as the TDE does, and no fallback to another
  estimator. `confidence` is the smallest maximum of the reference's cross-correlations; `events` stays 0.
- **blobmodel ties the periodic cell to the grid**, so `make_2d_realization` sets `geometry.Ly` after
  `Geometry.from_arrays` (10 by default, at least the view plus three widths on each side, the mirror copies then
  changing the correlations by less than $e^{-9}$). This couples the helper to blobmodel's internals.
- **Figures** in fusion_scripts go through `apd_check/figures.save` (figure_provenance, `stacklevel=2`), from
  `@provenance` functions with current docstrings, regenerated by the module's entry point. The field plots colour the
  confidence from Okabe–Ito blue to vermillion over the deck's range (0.5 to 1, `extend="min"`), use one arrow scale
  (2100 m/s per cm) on every shot, grey crosses for dead pixels and grey circles for live pixels without an estimate.
  `import netCDF4` before anything that imports h5py.
- **Provenance of the FCA figures**: `commit` 3f97362 with a non-empty `uncommitted` list, the author's
  letter-presentation files that were in progress at the time (committed since in fusion_scripts `d02e850`;
  `config.py` gained only a path for them). None of them affects the fields; regenerating would record a clean tree.
- **fusion_scripts' venv** holds the branch's 1.1.0 (`64260fd`) in place of PyPI's 1.0.0. The clean 1.1.0 wheel has no
  `conditional_averaging` module (dropped from the repository in 2023; only a stale `build/` directory kept it in the
  1.0.0 wheel), which `decorrelation/apd_check/simulation_column.py` imports.
- **Paper text** (the note, the letter) follows the author's rules in `~/Dropbox/Obsidian Vault/Papers/paper_style.md`.

## Literature on the conditioning (2026-11-10)

Compiled for the note's Limitations bullet; the bibliographic data of every entry verified through Crossref.
*Read* means the text was read; *abstract*, only the abstract (publisher pages refuse automated fetches, so abstracts
came from Crossref and OpenAlex).

- **Briggs 1984** (MAP Handbook 13, *read*; `briggs_analysis_1984`): place the three antennas so that the zero-lag
  cross-correlation is of the order of 0.5; closer, the "true" velocity comes out too small, and farther apart the
  analysis frequently breaks down; it is "undesirable to work with very high values of correlation" (no derivation;
  the $(\ell/\Delta)^2$ amplification is our reading). Unweighted least squares over all pairs; a fit to the complete
  cross-correlation functions is called preferable but not in routine use. The crossing times are "much less"
  uncertain than the lags of the maxima. The zero-lag spike is removed by interpolating the autocorrelation through the
  points on each side and renormalizing: the spaced-antenna version of the package's spike-free normalization. Rejects
  hyperbolic contours ($W^{-1}$ not positive definite) and apparent and true velocities that differ greatly.
- **Triangle size effect**, the FCA speed falling with decreasing spacing: Meek 1990 (*abstract*) takes noise as the
  only significant factor and removes the discrepancy by removing the noise term; Holdsworth and Reid 1995
  (*abstract*) confirm that the major cause is the failure to compensate for noise; Holdsworth 1999 (*abstract*) adds
  coarse digitization, receiver differences and antenna coupling, and further instrumental effects on simulated data.
  Golley and Rossiter 1970 observed it first (*secondhand*, not cited). Burke 1987 (*abstract*) gives a model cause,
  random surfaces that change along their motion, which is not the pulse model; not cited.
- **Fits to all lags**: Fedor 1967 and Baker 1981 (*abstracts*), least squares to the correlation functions of any
  number of sensors, the literature's way of using more of the data. No paper found that weights the zero-lag
  equations by $R$ or its variance.
- **Plasma imaging**: Zweben et al. 2017 (*read*) fit the NSTX gas-puff-imaging zero-lag correlation to an ellipse at
  the contour $C = 0.8$, several pixels from the origin rather than at the nearest neighbours, and reject axes below
  one pixel: the analogue of a larger `max_separation`, a different estimator; not cited. Guszejnov et al. 2013
  estimate the tilt from time delays; nothing on zero-lag conditioning. No gas-puff-imaging, beam-emission or probe
  paper found on the conditioning of a zero-lag ellipse fit.
- Not reached: May 1988 and Hocking et al. 1989 (no abstract available), Holdsworth's thesis, the 1976 Canadian
  report on Golley and Rossiter.

## What is left

1. Push `fca`, merge into `main` and release 1.1.0 (PyPI has 1.0.0; readthedocs reads `docs/requirements.txt`, which
   now asks for `blobmodel>=2.1.0`).
2. On C-Mod: `max_separation = 2` for the conditioning of the shape fit, above all below $f_\mathrm{GW} = 0.65$; the
   2DCA comparison at low density; the sign flip of the poloidal arrows next to the separatrix.
3. Only if data show the need, the simplest version first: weighting the shape equations by $R$, Briggs's calibration
   against the autocorrelation, G-INT, the lag-dependent shape with a velocity spread (out of scope in the plan).
4. The slides (`30c2413`) predate the note's conditioning bullet (`b612cc7`): add it to their limitations when next
   rendered.
5. fusion_scripts: `simulation_column.py` and the `conditional_averaging` gap; `fca_fields.py`'s module docstring says
   the FCA gave "half" the 3TDE speed on shot 027, where the measured ratio is 0.63; regenerating the figures clears the
   `uncommitted` list of their provenance.
6. `test_cond_av` stays marginal, as before the branch.
