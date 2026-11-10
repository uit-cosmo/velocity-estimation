# Slides of the FCA note

Manim slides of the note in `../note/main.tex`: velocity estimation by the crossing time and full correlation analysis
(FCA), as implemented in `velocity_estimation`. They cover:

- a superposition of pulses and its cross-correlation function;
- the crossing time, why the decorrelation cancels, and the crossing equation;
- the shape $W^{-1}$ from the zero-lag correlations, and the FCA velocity;
- the pooled least-squares fit around each reference pixel;
- the correlations from measured signals: the noise spike and the sampled crossing;
- the usage and the limitations.

Every figure is drawn in manim; there are no assets to prepare.

## Install, once

On Linux, manim builds against the cairo and pango headers. Manim 0.19.1 is the last release that supports Python 3.10.

```bash
sudo apt install libcairo2-dev libpango1.0-dev pkg-config
python3 -m venv .venv
.venv/bin/pip install "manim==0.19.1" "manim-slides[pyside6]"
```

The venv of `fusion_scripts/decorrelation/letter_presentation` has the same versions and works as well.

## Render and present

Run these from this folder. The scenes are listed in `SCENES` in `slides.py`. Use `--quality=l` for a quick draft. Write
the quality out in full: `manim-slides render` reads the `h` of `-qh` as its own help flag. Render the scenes one after
the other, not in parallel processes, which race on manim's LaTeX cache.

```bash
S="Title Superposition Correlation CrossingTime CrossingEquation Shape FCA Grid FromData Usage Summary"
.venv/bin/manim-slides render --quality=h slides.py $S
.venv/bin/manim-slides present $S
.venv/bin/manim-slides convert $S talk.html
```

The last command writes a reveal.js page that you can open in a browser.
