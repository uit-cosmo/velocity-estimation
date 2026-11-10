"""Manim slides of the note: velocity estimation by the crossing time and FCA.

The talk follows ../note/main.tex: from a superposition of pulses to the
cross-correlation, the crossing time, the shape from the zero-lag correlations,
FCA, the pooled fit on the imaging grid, the correlations from measured signals,
the usage and the limitations. README.md has the commands to render and present.
"""

import numpy as np
from manim import *
from manim_slides import Slide

SCENES = ["Title", "Superposition", "Correlation", "CrossingTime", "CrossingEquation", "Shape", "FCA", "Grid",
          "FromData", "Usage", "Summary"]

config.background_color = WHITE

# Okabe-Ito colours
OI_BLUE, OI_ORANGE, OI_GREEN, OI_RED, OI_SKY, OI_PURPLE = (
    "#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7")
INK, MUTED = "#000000", "#777777"

# a wide line, so that a Tex string breaks only where it says \\
TEMPLATE = TexTemplate()
TEMPLATE.add_to_preamble(r"\usepackage{mathptmx}\usepackage{bm}"
                         r"\setlength{\textwidth}{40cm}\setlength{\linewidth}{40cm}")
Tex.set_default(color=INK, tex_template=TEMPLATE)
MathTex.set_default(color=INK, tex_template=TEMPLATE)
Dot.set_default(color=INK)
Line.set_default(color=INK)

AXES = {"color": INK, "stroke_width": 2, "include_tip": False,
        "decimal_number_config": {"num_decimal_places": 0, "color": INK}}


# ---------------------------------------------------------------- helpers

def lines(*rows, size=28, color=INK, buff=0.14, align=LEFT):
    return VGroup(*[Tex(r, font_size=size, color=color) for r in rows]).arrange(DOWN, aligned_edge=align, buff=buff)


def bullets(*items, size=28, buff=0.32):
    """Items are strings or tuples of lines, each with a dot before its first line."""
    rows = VGroup()
    for item in items:
        body = lines(*(item if isinstance(item, tuple) else (item,)), size=size)
        rows.add(VGroup(Dot(radius=0.045, color=OI_BLUE).next_to(body[0], LEFT, buff=0.22), body))
    return rows.arrange(DOWN, aligned_edge=LEFT, buff=buff)


def blob(width, height, color=OI_BLUE, angle=0.0, layers=7, opacity=0.16):
    """Concentric ellipses whose stacked opacity falls off like a Gaussian pulse."""
    ellipses = VGroup(*[Ellipse(width=width * s, height=height * s, stroke_width=0, fill_color=color,
                                fill_opacity=opacity) for s in np.linspace(1, 0.15, layers)])
    return ellipses.rotate(angle)


def fade(pulse, amplitude, opacity=0.16):
    for ellipse in pulse:
        ellipse.set_fill(opacity=opacity * amplitude)


def probe(point, label=None, size=26, side=DOWN):
    mark = Square(side_length=0.14, fill_color=INK, fill_opacity=1, stroke_width=0).move_to(point)
    if label is None:
        return mark
    return VGroup(mark, MathTex(label, font_size=size).next_to(mark, side, buff=0.12))


def w_inverse(l1, l2, alpha):
    """The inverse of W = Q diag(l1^2, l2^2) Q^T, the shape of the correlation."""
    q = np.array([[np.cos(alpha), -np.sin(alpha)], [np.sin(alpha), np.cos(alpha)]])
    return q @ np.diag([1 / l1 ** 2, 1 / l2 ** 2]) @ q.T


def contour(centre, through, l1, l2, alpha, color=OI_BLUE, **kwargs):
    """The contour of the correlation centred on `centre` that passes through the point `through`."""
    d = np.asarray(through[:2]) - np.asarray(centre[:2])
    k = max(d @ w_inverse(l1, l2, alpha) @ d, 1e-6)
    return Ellipse(width=2 * l1 * np.sqrt(k), height=2 * l2 * np.sqrt(k), color=color, **kwargs
                   ).rotate(alpha).move_to([centre[0], centre[1], 0])


class Talk(Slide):
    """A scene whose slides each end on the final state of their last animation."""

    # reversed videos only serve to play a slide backwards, and their split
    # into segments hangs manim-slides 5.7 on long 1080p slides
    skip_reversing = True

    def construct(self):
        self.slides()
        self.wait(0.3)

    def next_slide(self, *args, **kwargs):
        # play() stops one frame short of the end, which a slide would freeze on
        self.wait(0.3)
        super().next_slide(*args, **kwargs)

    def heading(self, text):
        title = Tex(text, font_size=40).to_corner(UL, buff=0.45)
        y = title.get_bottom()[1] - 0.12
        half = config.frame_width / 2 - 0.45
        return VGroup(title, Line([-half, y, 0], [half, y, 0], color=OI_BLUE, stroke_width=2))

    def cite(self, text):
        return Tex(text, font_size=20, color=MUTED).to_corner(DR, buff=0.2)

    def clear(self):
        for mob in self.mobjects:
            mob.clear_updaters()
        if self.mobjects:
            self.play(*[FadeOut(mob) for mob in self.mobjects])


# ---------------------------------------------------------------- scenes

class Title(Talk):
    def slides(self):
        title = Tex(r"Velocity estimation by the crossing time", font_size=50)
        title2 = Tex(r"and full correlation analysis", font_size=50)
        sub = Tex(r"the model, the method and its implementation in \texttt{velocity\_estimation}",
                  font_size=30, color=MUTED)
        authors = Tex(r"J.~M.~Losada and O.~E.~Garcia", font_size=32)
        place = Tex(r"Department of Physics and Technology, UiT The Arctic University of Norway",
                    font_size=24, color=MUTED)
        VGroup(title, title2, sub, authors, place).arrange(DOWN, buff=0.3).shift(UP * 1.3)
        points = [np.array(p) for p in [(-0.8, -2.9, 0), (0.8, -2.9, 0), (-0.8, -1.5, 0)]]
        probes = VGroup(*[probe(p) for p in points])
        t = ValueTracker(0.0)
        pulse = blob(2.4, 0.9, angle=120 * DEGREES).move_to([-8.6, -2.2, 0])

        def follow(m):
            tt = t.get_value()
            m.move_to([-8.6 + 17.2 * tt, -2.2, 0])
            fade(m, np.exp(-((tt - 0.5) / 0.3) ** 2))

        pulse.add_updater(follow)
        self.play(Write(title), Write(title2), FadeIn(sub, shift=UP * 0.2))
        self.play(FadeIn(authors), FadeIn(place), FadeIn(probes))
        self.add(pulse)
        self.play(t.animate.set_value(1.0), run_time=4, rate_func=linear)


class Superposition(Talk):
    """Pulses that move with one velocity, fade in and out, and are sampled on a grid."""

    def slides(self):
        head = self.heading(r"A superposition of pulses")
        box = Rectangle(width=5.6, height=4.4, color=MUTED, stroke_width=1.5).move_to([-3.6, -0.55, 0])
        x0, y0 = box.get_left()[0], box.get_bottom()[1]
        grid = [np.array([x0 + 1.0 + 1.2 * i, y0 + 0.8 + 0.95 * j, 0]) for i in range(4) for j in range(4)]
        probes = VGroup(*[probe(p) for p in grid])
        watched = grid[6]
        ring = Circle(radius=0.2, color=OI_RED, stroke_width=3).move_to(watched)

        rng = np.random.default_rng(3)
        v = np.array([1.0, -0.35])
        l1, l2, alpha, tau_p, duration = 0.55, 0.25, 30 * DEGREES, 1.1, 8.0
        w_inv = w_inverse(l1, l2, alpha)
        n = 40
        starts = rng.uniform(-1.0, duration + 1.0, n)
        centres = np.column_stack([rng.uniform(x0 - 0.5, x0 + 5.6, n), rng.uniform(y0 + 0.3, y0 + 4.6, n)])
        amps = rng.exponential(1.0, n).clip(0.3, 2.0)
        t = ValueTracker(0.0)

        def position(k, tt):
            return centres[k] + v * (tt - starts[k])

        def envelope(k, tt):
            return amps[k] * np.exp(-((tt - starts[k]) / tau_p) ** 2)

        def inside(p):
            # fade the pulses out at the edge of the view
            dx = min(p[0] - box.get_left()[0], box.get_right()[0] - p[0])
            dy = min(p[1] - box.get_bottom()[1], box.get_top()[1] - p[1])
            return float(np.clip(min(dx, dy) / 0.4, 0, 1))

        pulses = VGroup()
        for k in range(n):
            pulse = blob(4 * l1, 4 * l2, angle=alpha)

            def follow(m, k=k):
                p = position(k, t.get_value())
                m.move_to([p[0], p[1], 0])
                fade(m, min(1.0, envelope(k, t.get_value()) / 1.5) * inside(p), opacity=0.2)

            pulse.add_updater(follow)
            pulses.add(pulse)

        def signal(tt):
            total = 0.0
            for k in range(n):
                d = watched[:2] - position(k, tt)
                total += envelope(k, tt) * np.exp(-d @ w_inv @ d)
            return total

        axes = Axes(x_range=[0, duration, 2], y_range=[0, 2.2, 1], x_length=5.6, y_length=2.0,
                    axis_config=AXES).move_to([3.4, -1.9, 0])
        labels = VGroup(MathTex(r"t", font_size=28).next_to(axes.x_axis, RIGHT, buff=0.1),
                        MathTex(r"\Phi(\bm x_\text{p}, t)", font_size=26, color=OI_RED
                                ).next_to(axes.y_axis, UP, buff=0.1))
        trace = always_redraw(lambda: axes.plot(signal, x_range=[0, max(0.02, t.get_value()), 0.02],
                                                color=OI_RED, stroke_width=2.5))
        equation = MathTex(r"\Phi(\bm x,t) = \sum_{k=1}^{K} a_k\, \psi\big(\bm x - \bm\xi_k - \bm v(t-s_k)\big)\,"
                           r"\varphi\!\left(\frac{t-s_k}{\tau_\text{p}}\right)", font_size=32).move_to([3.4, 1.6, 0])
        legend = bullets(r"$K$ uncorrelated pulses, amplitudes $a_k$",
                         r"all move with the same velocity $\bm v$",
                         r"shape $\psi$, envelope $\varphi$ with lifetime $\tau_\text{p}$",
                         size=24, buff=0.18).next_to(equation, DOWN, buff=0.35)
        self.play(FadeIn(head), Create(box), FadeIn(probes))
        self.add(pulses)
        self.play(Write(equation))
        self.play(Create(ring), FadeIn(axes), FadeIn(labels))
        self.add(trace)
        self.play(t.animate.set_value(duration), run_time=10, rate_func=linear)
        self.next_slide()
        self.play(FadeIn(legend, lag_ratio=0.3))
        self.next_slide()
        self.clear()

        # the pulse shape
        head = self.heading(r"Elongated, tilted pulses")
        theta = 30 * DEGREES
        centre = np.array([-3.6, -0.3, 0])
        ellipse = Ellipse(width=4.4, height=1.6, color=OI_BLUE, fill_opacity=0.15).rotate(theta).move_to(centre)
        major = np.array([np.cos(theta), np.sin(theta), 0])
        minor = np.array([-np.sin(theta), np.cos(theta), 0])
        axis1 = Line(centre, centre + 2.2 * major, stroke_width=3)
        axis2 = Line(centre, centre + 0.8 * minor, stroke_width=3)
        horizontal = DashedLine(centre + 2.8 * LEFT, centre + 2.8 * RIGHT, color=MUTED)
        angle = Angle(horizontal, axis1, radius=1.0, other_angle=False)
        names = VGroup(MathTex(r"\ell_1", font_size=32).next_to(axis1.get_center(), UL, buff=0.08),
                       MathTex(r"\ell_2", font_size=32).next_to(axis2.get_center(), LEFT, buff=0.1),
                       MathTex(r"\alpha", font_size=30).move_to(centre + 1.3 * np.array([np.cos(theta / 2),
                                                                                       np.sin(theta / 2), 0])))
        velocity = Arrow(centre, centre + np.array([2.2, -1.2, 0]), buff=0, color=OI_ORANGE, stroke_width=5)
        vname = MathTex(r"\bm v", font_size=32, color=OI_ORANGE).next_to(velocity.get_end(), DR, buff=0.05)
        equations = VGroup(
            MathTex(r"\psi(\bm x) \propto \exp\left(-\bm x^{\mathsf T} W^{-1}\bm x\right),\qquad"
                    r"\varphi(\theta)\propto e^{-\theta^2}", font_size=32),
            MathTex(r"W = Q_\alpha\begin{pmatrix}\ell_1^2 & 0\\ 0 & \ell_2^2\end{pmatrix}Q_\alpha^{\mathsf T}",
                    font_size=34),
            MathTex(r"w_{xy} = \left(\frac{1}{\ell_1^2}-\frac{1}{\ell_2^2}\right)\sin\alpha\cos\alpha", font_size=32),
            Tex(r"$w_{ij} = (W^{-1})_{ij}$; $w_{xy}=0$ for round or axis-aligned pulses", font_size=24, color=MUTED),
        ).arrange(DOWN, buff=0.35).move_to([3.3, 0.0, 0])
        self.play(FadeIn(head), DrawBorderThenFill(ellipse), Create(horizontal))
        self.play(Create(axis1), Create(axis2), Create(angle), FadeIn(names))
        self.play(GrowArrow(velocity), FadeIn(vname))
        self.play(Write(equations[0]))
        self.play(Write(equations[1]))
        self.play(FadeIn(equations[2:]))


class Correlation(Talk):
    """The cross-correlation: the shape centred on v tau, lowered by the decorrelation."""

    def slides(self):
        head = self.heading(r"The cross-correlation function")
        origin = np.array([-4.2, 0.2, 0])
        plane = VGroup(Arrow(origin + 2.4 * LEFT, origin + 3.4 * RIGHT, buff=0, color=MUTED, stroke_width=2,
                             tip_length=0.18),
                       Arrow(origin + 2.6 * DOWN, origin + 2.4 * UP, buff=0, color=MUTED, stroke_width=2,
                             tip_length=0.18))
        plane_labels = VGroup(MathTex(r"\Delta_x", font_size=28).next_to(plane[0], RIGHT, buff=0.08),
                              MathTex(r"\Delta_y", font_size=28).next_to(plane[1], UP, buff=0.08))
        v = np.array([1.0, -0.45, 0])
        tau_p = 1.3
        tau = ValueTracker(0.0)
        path = DashedLine(origin - 0.6 * v, origin + 3.2 * v, color=OI_ORANGE)
        shape = blob(2.6, 1.0, angle=30 * DEGREES, opacity=0.2)

        def follow(m):
            s = tau.get_value()
            m.move_to(origin + s * v)
            fade(m, np.exp(-s ** 2 / (2 * tau_p ** 2)), opacity=0.2)

        shape.add_updater(follow)
        centre = always_redraw(lambda: Dot(origin + tau.get_value() * v, radius=0.05, color=INK))
        readout = always_redraw(lambda: MathTex(r"\tau = %.1f" % tau.get_value(), font_size=28
                                                ).move_to(origin + np.array([-1.4, 2.0, 0])))
        vt = MathTex(r"\bm v\tau", font_size=28, color=OI_ORANGE).move_to(origin + 2.6 * v + 0.45 * DOWN)
        equations = VGroup(
            MathTex(r"R(\bm\Delta,\tau) = \frac{\langle \delta\Phi(\bm x,t)\,\delta\Phi(\bm x+\bm\Delta,t+\tau)\rangle}"
                    r"{\langle\delta\Phi^2\rangle}", font_size=32),
            MathTex(r"R(\bm\Delta,\tau) = \rho_\psi(\bm\Delta-\bm v\tau)\,", r"\rho_\varphi(\tau/\tau_\text{p})",
                    font_size=36),
        ).arrange(DOWN, buff=0.45).move_to([3.4, 1.4, 0])
        equations[1][1].set_color(OI_RED)
        campbell = Tex(r"Campbell's theorem", font_size=24, color=MUTED).next_to(equations[1], UP, buff=0.12)
        notes = bullets(r"at each lag: the pulse autocorrelation centred on $\bm v\tau$",
                        r"decorrelation: a factor of the lag alone",
                        size=26, buff=0.2).next_to(equations, DOWN, buff=0.5)
        notes[1][1].set_color(OI_RED)
        gauss = MathTex(r"R = e^{-\tau^2/2\tau_\text{p}^2}\,\exp\left[-\tfrac12(\bm\Delta-\bm v\tau)^{\mathsf T}W^{-1}"
                        r"(\bm\Delta-\bm v\tau)\right]", font_size=30, color=OI_BLUE).to_edge(DOWN, buff=0.45)
        gauss_note = Tex(r"Gaussian pulses: $W$ is the covariance of the correlation", font_size=24, color=MUTED
                         ).next_to(gauss, UP, buff=0.12)
        self.play(FadeIn(head), Create(plane), FadeIn(plane_labels))
        self.play(Write(equations[0]))
        self.play(FadeIn(campbell), Write(equations[1]))
        self.add(shape, centre, readout)
        self.play(Create(path), FadeIn(vt))
        self.play(tau.animate.set_value(2.6), run_time=5, rate_func=linear)
        self.play(FadeIn(notes, lag_ratio=0.3))
        self.next_slide()
        self.play(FadeIn(gauss_note), Write(gauss))


class CrossingTime(Talk):
    """The autocorrelation equals the cross-correlation at the crossing time, whatever the lifetime."""

    def slides(self):
        head = self.heading(r"The crossing time")
        axes = Axes(x_range=[-0.5, 2, 0.5], y_range=[0, 1.05, 0.5], x_length=6.4, y_length=3.8,
                    axis_config=AXES, x_axis_config={"numbers_to_include": [1, 1.5, 2],
                                                     "decimal_number_config": {"num_decimal_places": 1,
                                                                               "color": INK},
                                                     "font_size": 22}
                    ).move_to([-2.8, -0.6, 0])
        tau_label = MathTex(r"\tau v/\Delta", font_size=28).next_to(axes.x_axis, DOWN, buff=0.45)
        tau_p = ValueTracker(20.0)

        def auto(s):
            return np.exp(-s ** 2 / (2 * tau_p.get_value() ** 2) - s ** 2 / 2)

        def cross(s):
            return np.exp(-s ** 2 / (2 * tau_p.get_value() ** 2) - (1 - s) ** 2 / 2)

        curves = always_redraw(lambda: VGroup(
            axes.plot(auto, x_range=[-0.5, 2], color=INK, stroke_width=3),
            DashedVMobject(axes.plot(cross, x_range=[-0.5, 2], color=OI_BLUE, stroke_width=3), num_dashes=40)))
        meet = always_redraw(lambda: Dot(axes.c2p(0.5, auto(0.5)), radius=0.09, color=OI_GREEN))
        drop = always_redraw(lambda: DashedLine(axes.c2p(0.5, 0), axes.c2p(0.5, auto(0.5)), color=OI_GREEN))
        names = VGroup(MathTex(r"R(\bm 0,\tau)", font_size=28).move_to(axes.c2p(-0.15, 1.12)),
                       MathTex(r"R(\bm\Delta,\tau)", font_size=28, color=OI_BLUE).move_to(axes.c2p(1.55, 1.12)),
                       MathTex(r"\tau_\text{e}", font_size=28, color=OI_GREEN).next_to(axes.c2p(0.5, 0), DOWN,
                                                                                         buff=0.12))
        definition = VGroup(Tex(r"$\tau_\text{e}$: the lag at which", font_size=30),
                            MathTex(r"R(\bm 0,\tau_\text{e}) = R(\bm\Delta,\tau_\text{e})", font_size=36)
                            ).arrange(DOWN, buff=0.18).move_to([3.7, 2.0, 0])
        lifetime = DecimalNumber(20.0, num_decimal_places=2, font_size=32, color=MUTED)
        lifetime.add_updater(lambda m: m.set_value(tau_p.get_value()))
        knob = VGroup(MathTex(r"\tau_\text{p} v/\Delta =", font_size=32, color=MUTED), lifetime
                      ).arrange(RIGHT, buff=0.15).move_to([3.7, 0.6, 0])
        ve = VGroup(MathTex(r"\widehat v_\text{e} = \frac{\Delta}{2\tau_\text{e}} = v", font_size=38,
                            color=OI_GREEN)).move_to([3.7, -0.5, 0])
        why = lines(r"both sides carry the same $\rho_\varphi(\tau/\tau_\text{p})$:",
                    r"the decorrelation cancels", size=28, color=OI_BLUE).move_to([3.7, -1.9, 0])
        self.play(FadeIn(head), FadeIn(axes), FadeIn(tau_label), FadeIn(definition))
        self.add(curves)
        self.play(FadeIn(names[:2]))
        self.add(drop, meet)
        self.play(FadeIn(names[2]), Write(ve))
        self.next_slide()
        self.play(FadeIn(knob))
        self.play(tau_p.animate.set_value(0.45), run_time=5)
        self.play(FadeIn(why))
        self.next_slide()
        self.clear()

        # the geometry for round pulses
        head = self.heading(r"Where it crosses: round pulses")
        unit, beta = 2.0, 40 * DEGREES
        p0 = np.array([-4.8, -1.6, 0])
        p1 = p0 + 2 * unit * RIGHT
        e_v = np.array([np.cos(beta), np.sin(beta), 0])
        probes = VGroup(probe(p0, r"P_0", side=DL), probe(p1, r"P_1", side=DR))
        separation = Line(p0, p1, stroke_width=2)
        bisector = DashedLine(p0 + unit * RIGHT + 0.6 * DOWN, p0 + unit * RIGHT + 3.6 * UP, color=MUTED)
        path = DashedLine(p0 - 0.3 * e_v, p0 + 3.6 * e_v, color=OI_ORANGE)
        tau = ValueTracker(0.05)
        tau_e = unit / np.cos(beta)

        def centre():
            return p0 + tau.get_value() * e_v

        circle = always_redraw(lambda: Circle(radius=np.linalg.norm(centre() - p0), color=OI_BLUE,
                                              stroke_width=3).move_to(centre()))
        arrow = always_redraw(lambda: Arrow(p0, centre(), buff=0, color=OI_ORANGE, stroke_width=5,
                                            max_tip_length_to_length_ratio=0.15))
        vt = always_redraw(lambda: MathTex(r"\bm v\tau", font_size=30, color=OI_ORANGE
                                           ).next_to(centre(), UL, buff=0.06))
        text = VGroup(
            Tex(r"contour of $R(\cdot,\tau)$ through $P_0$:", font_size=28),
            Tex(r"it reaches $P_1$ at $\tau = \tau_\text{e}$", font_size=28),
            MathTex(r"|\bm v\tau_\text{e}| = |\bm\Delta - \bm v\tau_\text{e}|", font_size=34),
            MathTex(r"\Rightarrow\quad \hat{\bm e}\cdot\bm v = \frac{\Delta}{2\tau_\text{e}}", font_size=38,
                    color=OI_GREEN),
            Tex(r"the component along the separation,", font_size=26, color=MUTED),
            Tex(r"for any direction of motion and any lifetime", font_size=26, color=MUTED),
        ).arrange(DOWN, buff=0.22).move_to([3.4, 0.2, 0])
        self.play(FadeIn(head), FadeIn(probes), Create(separation), Create(bisector), Create(path))
        self.add(circle, arrow, vt)
        self.play(FadeIn(text[:2]))
        self.play(tau.animate.set_value(tau_e), run_time=4, rate_func=linear)
        self.play(Flash(p1, color=OI_GREEN), Indicate(probes[1], color=OI_GREEN))
        self.play(Write(text[2]))
        self.play(Write(text[3]), FadeIn(text[4:]))


class CrossingEquation(Talk):
    """The crossing equation in general, and what it gives for elongated pulses."""

    def slides(self):
        head = self.heading(r"The crossing equation")
        steps = VGroup(
            Tex(r"Elliptic shape: $\rho_\psi(\bm r) = g(\bm r^{\mathsf T}W^{-1}\bm r)$, $g$ decreasing", font_size=30),
            MathTex(r"g\big((\bm v\tau_\text{e})^{\mathsf T}W^{-1}\bm v\tau_\text{e}\big) = "
                    r"g\big((\bm\Delta-\bm v\tau_\text{e})^{\mathsf T}W^{-1}(\bm\Delta-\bm v\tau_\text{e})\big)",
                    font_size=32),
            MathTex(r"\bm\Delta^{\mathsf T}W^{-1}\bm v = \frac{\bm\Delta^{\mathsf T}W^{-1}\bm\Delta}{2\tau_\text{e}}",
                    font_size=42, color=OI_GREEN),
            Tex(r"linear in $\bm v$; any envelope, any elliptic shape;"
                r" no crossing $\Leftrightarrow$ $1/\tau_\text{e}=0$", font_size=28, color=MUTED),
        ).arrange(DOWN, buff=0.4).next_to(head, DOWN, buff=0.55)
        box = SurroundingRectangle(steps[2], color=OI_GREEN, buff=0.15)
        round_case = MathTex(r"W=\ell^2 I:\quad \frac{\Delta}{2\tau_\text{e}} = \hat{\bm e}\cdot\bm v", font_size=34
                             ).next_to(steps, DOWN, buff=0.5)
        tilted_case = MathTex(r"\text{otherwise:}\quad \frac{\Delta}{2\tau_\text{e}} = "
                              r"\frac{\hat{\bm e}^{\mathsf T}W^{-1}\bm v}{\hat{\bm e}^{\mathsf T}W^{-1}\hat{\bm e}}"
                              r"\neq \hat{\bm e}\cdot\bm v", font_size=34, color=OI_RED
                              ).next_to(round_case, DOWN, buff=0.3)
        self.play(FadeIn(head), FadeIn(steps[0]))
        self.play(Write(steps[1]))
        self.play(Write(steps[2]), Create(box))
        self.play(FadeIn(steps[3]))
        self.next_slide()
        self.play(Write(round_case))
        self.play(Write(tilted_case))
        self.next_slide()
        self.clear()

        # elongated pulses: an apparent velocity across the motion
        head = self.heading(r"Elongated, tilted pulses moving along $x$")
        unit = 1.25
        l1, l2, alpha, delta = 1.4 * unit, 0.55 * unit, 120 * DEGREES, 1.6 * unit
        p0 = np.array([-4.6, -0.9, 0])
        p1 = p0 + delta * UP
        w_inv = w_inverse(l1, l2, alpha)
        s_e = delta * w_inv[1, 1] / (2 * w_inv[0, 1])
        probes = VGroup(probe(p0, r"P_0", side=DL), probe(p1, r"P_1", side=LEFT))
        path = DashedLine(p0 + 1.4 * LEFT, p0 + 3.4 * RIGHT, color=OI_ORANGE)
        bisector = DashedLine(p0 + delta / 2 * UP + 1.4 * LEFT, p0 + delta / 2 * UP + 3.4 * RIGHT, color=MUTED,
                              dash_length=0.05)
        normal = w_inv @ np.array([0, delta])
        direction = np.array([-normal[1], normal[0], 0]) / np.linalg.norm(normal)
        w_bisector = DashedLine(p0 + delta / 2 * UP - 2.4 * direction, p0 + delta / 2 * UP + 2.4 * direction,
                                color=OI_PURPLE)
        names = VGroup(Tex(r"perpendicular bisector", font_size=22, color=MUTED).next_to(bisector, RIGHT, buff=0.1),
                       Tex(r"bisector in the metric of $W^{-1}$", font_size=22, color=OI_PURPLE
                           ).move_to([-6.6, -2.9, 0], aligned_edge=LEFT))
        tau = ValueTracker(0.05)

        def centre():
            return p0 + tau.get_value() * RIGHT

        ellipse = always_redraw(lambda: contour(centre() - p0, [0, 0], l1, l2, alpha, stroke_width=3
                                                ).shift(p0))
        arrow = always_redraw(lambda: Arrow(p0, centre(), buff=0, color=OI_ORANGE, stroke_width=5,
                                            max_tip_length_to_length_ratio=0.25))
        text = VGroup(
            MathTex(r"\bm v = v\hat{\bm x}", font_size=32, color=OI_ORANGE),
            Tex(r"the contour through $P_0$ still reaches $P_1$:", font_size=28),
            MathTex(r"\frac{\Delta}{2\tau_\text{e}} = v\,\frac{w_{xy}}{w_{yy}} \neq 0 = v_y", font_size=38,
                    color=OI_RED),
            Tex(r"an apparent velocity across the motion", font_size=28, color=OI_RED),
            Tex(r"$\Rightarrow$ measure $W^{-1}$: FCA", font_size=32, color=OI_GREEN),
        ).arrange(DOWN, buff=0.3).move_to([3.6, 0.0, 0])
        self.play(FadeIn(head), FadeIn(probes), Create(path), Create(bisector), FadeIn(names[0]))
        self.play(FadeIn(text[:2]))
        self.add(ellipse, arrow)
        self.play(tau.animate.set_value(s_e), run_time=4, rate_func=linear)
        self.play(Flash(p1, color=OI_RED), Create(w_bisector), FadeIn(names[1]))
        self.play(Write(text[2]), FadeIn(text[3]))
        self.play(FadeIn(text[4]))


class Shape(Talk):
    """The shape from the zero-lag correlations of three points."""

    def slides(self):
        head = self.heading(r"The shape from the zero-lag correlations")
        unit = 1.5
        p0 = np.array([-5.6, -2.2, 0])
        px, py = p0 + unit * RIGHT, p0 + unit * UP
        probes = VGroup(probe(p0, r"P_0", side=DL), probe(px, r"P_x", side=DR), probe(py, r"P_y", side=UL))
        pairs = VGroup(Line(p0, px, color=OI_ORANGE, stroke_width=4), Line(p0, py, color=OI_ORANGE, stroke_width=4),
                       Line(px, py, color=OI_ORANGE, stroke_width=4))
        zero = MathTex(r"-2\ln R(\bm\Delta,0) = \bm\Delta^{\mathsf T}W^{-1}\bm\Delta"
                       r" = \Delta_x^2 w_{xx} + 2\Delta_x\Delta_y w_{xy} + \Delta_y^2 w_{yy}", font_size=32
                       ).next_to(head, DOWN, buff=0.4)
        note = Tex(r"at zero lag the decorrelation is $\rho_\varphi(0)=1$; Gaussian pulses", font_size=24,
                   color=MUTED).next_to(zero, DOWN, buff=0.12)
        origin = np.array([0.0, -1.2, 0])
        scale = 1.25
        plane = VGroup(Line(origin + 2.6 * LEFT, origin + 2.6 * RIGHT, color=MUTED, stroke_width=2),
                       Line(origin + 2.4 * DOWN, origin + 2.4 * UP, color=MUTED, stroke_width=2))
        level = Ellipse(width=2 * 1.4 * scale, height=2 * 0.55 * scale, color=OI_BLUE, fill_opacity=0.12
                        ).rotate(120 * DEGREES).move_to(origin)
        w_inv = w_inverse(1.4, 0.55, 120 * DEGREES)
        dots, rays = VGroup(), VGroup()
        for d in ([1, 0], [0, 1], [-1, 1]):
            e = np.array(d, float) / np.linalg.norm(d)
            r = scale / np.sqrt(e @ w_inv @ e)
            point = origin + r * np.array([e[0], e[1], 0])
            rays.add(DashedLine(origin - 2.2 * np.array([e[0], e[1], 0]), origin + 2.2 * np.array([e[0], e[1], 0]),
                                color=MUTED))
            dots.add(VGroup(Dot(point, radius=0.08, color=OI_ORANGE), Dot(2 * origin - point, radius=0.08,
                                                                           color=OI_ORANGE)))
        level_name = MathTex(r"\bm\Delta^{\mathsf T}W^{-1}\bm\Delta=\kappa", font_size=28, color=OI_BLUE
                             ).move_to(origin + np.array([1.7, 1.6, 0]))
        triangle = VGroup(
            MathTex(r"w_{xx} = -\frac{2\ln R(\Delta\hat{\bm x},0)}{\Delta^2}", font_size=30),
            MathTex(r"w_{yy} = -\frac{2\ln R(\Delta\hat{\bm y},0)}{\Delta^2}", font_size=30),
            MathTex(r"w_{xy} = \frac{\ln R(\Delta(\hat{\bm y}-\hat{\bm x}),0) - \ln R(\Delta\hat{\bm x},0)"
                    r" - \ln R(\Delta\hat{\bm y},0)}{\Delta^2}", font_size=26),
        ).arrange(DOWN, buff=0.25, aligned_edge=LEFT).move_to([3.9, -0.6, 0])
        widths = Tex(r"eigenvalues $1/\ell_1^2,\,1/\ell_2^2$, eigenvector $\Rightarrow$ $\alpha$", font_size=26,
                     color=OI_GREEN).next_to(triangle, DOWN, buff=0.35)
        self.play(FadeIn(head), FadeIn(probes))
        self.play(Write(zero), FadeIn(note))
        self.play(Create(pairs))
        self.next_slide()
        self.play(Create(plane), FadeIn(rays))
        for k in range(3):
            self.play(Indicate(pairs[k], color=OI_RED), FadeIn(dots[k], scale=1.5), run_time=0.8)
        self.play(DrawBorderThenFill(level), FadeIn(level_name))
        self.next_slide()
        self.play(VGroup(plane, rays, dots, level, level_name).animate.scale(0.8).move_to([-1.7, -1.1, 0]))
        self.play(Write(triangle))
        self.play(FadeIn(widths))


class FCA(Talk):
    """The crossing equations with the measured shape, and the tilt scan."""

    def slides(self):
        head = self.heading(r"Full correlation analysis (FCA)")
        system = MathTex(r"\begin{pmatrix}\bm\Delta_1^{\mathsf T}W^{-1}\\ \bm\Delta_2^{\mathsf T}W^{-1}\end{pmatrix}"
                         r"\bm v = \begin{pmatrix}\bm\Delta_1^{\mathsf T}W^{-1}\bm\Delta_1/2\tau_{\text{e},1}\\"
                         r"\bm\Delta_2^{\mathsf T}W^{-1}\bm\Delta_2/2\tau_{\text{e},2}\end{pmatrix}", font_size=36
                         ).next_to(head, DOWN, buff=0.5)
        solution = MathTex(r"\widehat{\bm v}_\text{FCA} = W\begin{pmatrix}w_{xx}\,\Delta/2\tau_{\text{e},x}\\"
                           r"w_{yy}\,\Delta/2\tau_{\text{e},y}\end{pmatrix} = \bm v", font_size=40, color=OI_GREEN
                           ).next_to(system, DOWN, buff=0.45)
        remarks = bullets(r"crossing times: the decorrelation cancels",
                          r"zero-lag shape: back from the metric of $W^{-1}$ to that of the plane",
                          r"exact for Gaussian pulses of any elongation, tilt and lifetime, at any separation",
                          r"with $W^{-1}=I$: the crossing time of each pair (\texttt{CROSSING})",
                          size=27, buff=0.22).next_to(solution, DOWN, buff=0.5)
        self.play(FadeIn(head), Write(system))
        self.play(Write(solution))
        self.play(FadeIn(remarks, lag_ratio=0.3),
                  FadeIn(self.cite(r"Phillips and Spencer 1955; Briggs 1984; Holloway \textit{et al.} 1997")))
        self.next_slide()
        self.clear()

        # the tilt scan: CROSSING against FCA for pulses moving along x
        head = self.heading(r"Pulses moving along $x$, the tilt varied")
        alpha = ValueTracker(0.0)
        l1, l2, v = 1.0, 0.25, 1.0
        centre = np.array([-4.2, -0.3, 0])
        pulse = always_redraw(lambda: Ellipse(width=2.6, height=0.65, color=OI_BLUE, fill_opacity=0.15
                                              ).rotate(alpha.get_value()).move_to(centre))
        triangle = VGroup(probe(centre + np.array([-0.6, -2.4, 0])), probe(centre + np.array([0.4, -2.4, 0])),
                          probe(centre + np.array([-0.6, -1.4, 0])))
        origin = np.array([-0.5, -0.3, 0])
        scale = 1.6

        def crossing():
            w = w_inverse(l1, l2, alpha.get_value())
            return np.array([v, v * w[0, 1] / w[1, 1], 0])

        true = Arrow(origin, origin + scale * v * RIGHT, buff=0, color=INK, stroke_width=8)
        fca = DashedLine(origin, origin + scale * v * RIGHT, color=OI_GREEN, stroke_width=5)
        plain = always_redraw(lambda: Arrow(origin, origin + scale * crossing(), buff=0, color=OI_RED,
                                            stroke_width=5))
        legend = VGroup(
            VGroup(Line(ORIGIN, 0.6 * RIGHT, color=INK, stroke_width=8), Tex(r"$\bm v$", font_size=28)),
            VGroup(DashedLine(ORIGIN, 0.6 * RIGHT, color=OI_GREEN, stroke_width=5),
                   Tex(r"FCA", font_size=28, color=OI_GREEN)),
            VGroup(Line(ORIGIN, 0.6 * RIGHT, color=OI_RED, stroke_width=5),
                   Tex(r"crossing time, $W^{-1}=I$", font_size=28, color=OI_RED)),
        )
        for row in legend:
            row.arrange(RIGHT, buff=0.2)
        legend.arrange(DOWN, aligned_edge=LEFT, buff=0.2).move_to([4.3, 2.0, 0])
        angle = DecimalNumber(0, num_decimal_places=0, font_size=32, unit=r"^\circ", color=INK)
        angle.add_updater(lambda m: m.set_value(alpha.get_value() / DEGREES))
        readout = VGroup(MathTex(r"\alpha =", font_size=32), angle).arrange(RIGHT, buff=0.15
                                                                            ).next_to(pulse, UP, buff=1.0)
        note = Tex(r"aspect ratio $\ell_1/\ell_2 = 4$; the crossing time gives $v\,w_{xy}/w_{yy}$ along $y$",
                   font_size=26, color=MUTED).to_edge(DOWN, buff=0.45)
        self.play(FadeIn(head), FadeIn(pulse), FadeIn(triangle), FadeIn(readout))
        self.play(GrowArrow(true), Create(fca), FadeIn(legend))
        self.add(plain)
        self.play(FadeIn(note))
        self.play(alpha.animate.set_value(60 * DEGREES), run_time=4)
        self.play(alpha.animate.set_value(135 * DEGREES), run_time=5)


class Grid(Talk):
    """One pooled least-squares fit per reference pixel."""

    def slides(self):
        head = self.heading(r"On the imaging grid: one pooled fit per pixel")
        step = 1.1
        centre = np.array([-3.8, -0.6, 0])
        dots = VGroup(*[Dot(centre + step * np.array([i, j, 0]), radius=0.05, color=MUTED)
                        for i in range(-2, 3) for j in range(-2, 3)])
        p0 = centre
        neighbours = {"H_-": p0 + step * LEFT, "H_+": p0 + step * RIGHT, "V_+": p0 + step * UP,
                      "V_-": p0 + step * DOWN}
        sides = {"H_-": LEFT, "H_+": RIGHT, "V_+": UP, "V_-": DOWN}
        marks = VGroup(probe(p0, r"P_0", side=DR),
                       *[probe(p, name, side=sides[name]) for name, p in neighbours.items()])
        reference = VGroup(*[Line(p0, p, color=OI_ORANGE, stroke_width=6) for p in neighbours.values()])
        diagonals = VGroup(*[DashedLine(neighbours[h], neighbours[v], color=OI_BLUE, stroke_width=4)
                             for h in ("H_-", "H_+") for v in ("V_-", "V_+")])
        shape = VGroup(
            Tex(r"\textbf{Shape}: up to 8 equations, 3 unknowns", font_size=30),
            MathTex(r"\begin{pmatrix}\Delta_x^2 & 2\Delta_x\Delta_y & \Delta_y^2\end{pmatrix}_k"
                    r"\begin{pmatrix}w_{xx}\\w_{xy}\\w_{yy}\end{pmatrix} = -2\ln R(\bm\Delta_k,0)", font_size=30),
            Tex(r"$P_0$ with each neighbour, and each $H$ with each $V$", font_size=24, color=MUTED),
        ).arrange(DOWN, buff=0.2).move_to([3.3, 1.5, 0])
        velocity = VGroup(
            Tex(r"\textbf{Velocity}: up to 4 equations, 2 unknowns", font_size=30),
            MathTex(r"\bm\Delta_k^{\mathsf T}\widehat W^{-1}\bm v = "
                    r"\frac{\bm\Delta_k^{\mathsf T}\widehat W^{-1}\bm\Delta_k}{2\tau_{\text{e},k}}", font_size=32),
            Tex(r"$P_0$ with each neighbour; least squares", font_size=24, color=MUTED),
        ).arrange(DOWN, buff=0.2).move_to([3.3, -1.2, 0])
        notes = bullets(r"opposite neighbours keep the estimate centred on $P_0$",
                        r"a boundary pixel, one triangle: the same code, exactly determined",
                        r"$H_-H_+$ and $V_-V_+$ left out; only $H$ or only $V$: no estimate",
                        size=24, buff=0.14).to_edge(DOWN, buff=0.3).shift(RIGHT * 2.3)
        self.play(FadeIn(head), FadeIn(dots), FadeIn(marks))
        self.play(LaggedStart(*[Create(line) for line in reference], lag_ratio=0.25))
        self.play(LaggedStart(*[Create(line) for line in diagonals], lag_ratio=0.25))
        self.play(FadeIn(shape))
        self.next_slide()
        self.play(diagonals.animate.set_opacity(0.2), Indicate(reference, color=OI_RED))
        self.play(FadeIn(velocity))
        self.next_slide()
        self.play(FadeIn(notes, lag_ratio=0.3))


class FromData(Talk):
    """From measured signals to the inputs of the fit."""

    def slides(self):
        head = self.heading(r"Correlations from measured signals")
        dt = 0.2
        axes = Axes(x_range=[0, 1.05, 0.2], y_range=[0, 1.5, 0.5], x_length=5.4, y_length=3.4,
                    axis_config=AXES).move_to([-3.4, -0.4, 0])
        clean = axes.plot(lambda s: np.exp(-s ** 2 / 0.72), x_range=[0, 1.05], color=MUTED, stroke_width=2)
        samples = VGroup(*[Dot(axes.c2p(k * dt, np.exp(-(k * dt) ** 2 / 0.72)), radius=0.07)
                           for k in range(1, 6)])
        spike = Triangle(color=OI_RED, fill_opacity=1).scale(0.1).move_to(axes.c2p(0, 1.4))
        c1, c2 = np.exp(-dt ** 2 / 0.72), np.exp(-(2 * dt) ** 2 / 0.72)
        c0 = (4 * c1 - c2) / 3
        a = (c1 - c2) / (3 * dt ** 2)
        parabola = DashedVMobject(axes.plot(lambda s: c0 - a * s ** 2, x_range=[0, 2.2 * dt], color=OI_GREEN,
                                            stroke_width=4), num_dashes=12)
        star = Square(side_length=0.15, color=OI_GREEN, fill_opacity=1).rotate(45 * DEGREES
                                                                               ).move_to(axes.c2p(0, c0))
        labels = VGroup(MathTex(r"\tau", font_size=28).next_to(axes.x_axis, RIGHT, buff=0.1),
                        MathTex(r"C(\tau)", font_size=28).next_to(axes.y_axis, UP, buff=0.1),
                        Tex(r"white noise: a spike at zero lag", font_size=24, color=OI_RED
                            ).next_to(axes.c2p(0, 1.4), RIGHT, buff=0.25))
        formula = VGroup(
            MathTex(r"C^*(0) = \frac{4C(\mathrm{d}t) - C(2\,\mathrm{d}t)}{3}", font_size=36, color=OI_GREEN),
            MathTex(r"R_{ab}(\tau) = \frac{C_{ab}(\tau)}{\sqrt{C^*_{aa}(0)\,C^*_{bb}(0)}}", font_size=34),
            Tex(r"otherwise every $-2\ln R(\bm\Delta,0)$ gets a constant added,", font_size=24, color=MUTED),
            Tex(r"which biases $W^{-1}$; on by default (\texttt{remove\_noise\_spike})", font_size=24,
                color=MUTED),
        ).arrange(DOWN, buff=0.25).move_to([3.4, 0.4, 0])
        self.play(FadeIn(head), FadeIn(axes), FadeIn(labels[:2]), Create(clean))
        self.play(FadeIn(samples, lag_ratio=0.2))
        self.play(FadeIn(spike, scale=2), FadeIn(labels[2]))
        self.next_slide()
        self.play(Create(parabola), FadeIn(star))
        self.play(Write(formula[0]))
        self.play(Write(formula[1]), FadeIn(formula[2:]))
        self.next_slide()
        self.clear()

        # the sampled crossing
        head = self.heading(r"The sampled crossing")
        axes = Axes(x_range=[-1.2, 2.0, 0.5], y_range=[0, 1.05, 0.5], x_length=6.4, y_length=3.8,
                    axis_config=AXES).move_to([-3.0, -0.4, 0])

        def auto(s):
            return np.exp(-s ** 2 / 2.88 - s ** 2 / 2)

        def cross(s):
            return np.exp(-s ** 2 / 2.88 - (1 - s) ** 2 / 2)

        taus = np.arange(-1.2, 2.01, 0.2)
        auto_dots = VGroup(*[Circle(radius=0.06, color=INK, stroke_width=2).move_to(axes.c2p(s, auto(s)))
                             for s in taus])
        cross_dots = VGroup(*[Square(side_length=0.12, color=OI_BLUE, stroke_width=2).move_to(axes.c2p(s, cross(s)))
                              for s in taus])
        splines = VGroup(axes.plot(auto, x_range=[-1.2, 2.0], color=INK, stroke_width=1.5),
                         axes.plot(cross, x_range=[-1.2, 2.0], color=OI_BLUE, stroke_width=1.5))
        names = VGroup(MathTex(r"\bar R(\tau) = \tfrac12[R_{00}+R_{11}]", font_size=26).move_to(axes.c2p(-0.5, 1.12)),
                       MathTex(r"R(\bm\Delta,\tau)", font_size=26, color=OI_BLUE).move_to(axes.c2p(1.6, 0.95)))
        mark = DashedLine(axes.c2p(0.5, 0), axes.c2p(0.5, auto(0.5)), color=OI_GREEN)
        hit = Dot(axes.c2p(0.5, auto(0.5)), radius=0.09, color=OI_GREEN)
        te = MathTex(r"\tau_\text{e}", font_size=28, color=OI_GREEN).next_to(axes.c2p(0.5, 0), DOWN, buff=0.12)
        steps = bullets(r"first sign change of $\bar R - R(\bm\Delta,\cdot)$ on each side of zero lag",
                        r"refined on a cubic spline; the one nearest zero kept",
                        r"none within the window $[-T,T]$: $1/\tau_\text{e} = 0$, error $\le \Delta/2T$",
                        r"the mean autocorrelation: swapping the pixels flips the sign of $\tau_\text{e}$",
                        size=24, buff=0.22).move_to([3.5, -0.2, 0])
        self.play(FadeIn(head), FadeIn(axes))
        self.play(FadeIn(auto_dots, lag_ratio=0.05), FadeIn(cross_dots, lag_ratio=0.05), FadeIn(names))
        self.play(Create(splines))
        self.play(Create(mark), FadeIn(hit, scale=2), FadeIn(te))
        self.play(FadeIn(steps, lag_ratio=0.3))


class Usage(Talk):
    def slides(self):
        head = self.heading(r"Usage")
        source = "\n".join([
            "import velocity_estimation.two_dim_velocity_estimates as tdve",
            "from velocity_estimation import utils",
            "",
            "eo = tdve.EstimationOptions(velocity_method=tdve.VelocityMethod.FCA)",
            "eo.fca_options.window = None              # None: 100 dt",
            "eo.fca_options.remove_noise_spike = True",
            "eo.fca_options.minimum_correlation = 0.1  # shape fit",
            "",
            "md = tdve.estimate_velocity_field(utils.CModImagingDataInterface(ds), eo)",
            "vx, vy = md.get_vx(), md.get_vy()",
            "l1, l2, tilt = md.get_l1(), md.get_l2(), md.get_tilt()",
        ])
        code = Code(code_string=source, language="python", formatter_style="friendly", add_line_numbers=False,
                    background="rectangle", background_config={"fill_color": "#F4F4F4", "stroke_color": MUTED},
                    paragraph_config={"font": "Monospace", "font_size": 20})
        for part in code.family_members_with_points():
            if part is not code.background and color_to_rgb(part.get_fill_color()).min() > 0.95:
                part.set_fill(INK)
        code.next_to(head, DOWN, buff=0.45)
        notes = bullets(r"\texttt{VelocityMethod.CROSSING}: the crossing time pooled over the pairs; $\ell_1$, $\ell_2$,"
                        r" tilt are NaN",
                        r"the time-delay options are ignored; the default stays \texttt{VelocityMethod.TDE}",
                        r"the tilt is in radians from the first grid coordinate",
                        size=26, buff=0.2).next_to(code, DOWN, buff=0.45)
        self.play(FadeIn(head), FadeIn(code))
        self.play(FadeIn(notes, lag_ratio=0.3))


class Summary(Talk):
    def slides(self):
        head = self.heading(r"Summary")
        y0, dy = 2.0, 0.75
        xs = (-6.4, -4.4, -0.6, 1.9)
        titles = [Tex(t, font_size=26).move_to([x, y0, 0], aligned_edge=LEFT) for t, x in
                  zip((r"Method", r"Solves", r"Gives", r"Holds for"), xs)]
        table = VGroup(*titles, Line([-6.5, y0 - 0.3, 0], [6.7, y0 - 0.3, 0], stroke_width=1.5))
        rows = [(r"\texttt{CROSSING}", r"crossing equation, $W^{-1}=I$", r"$\hat{\bm e}\cdot\bm v$ per pair",
                 r"round pulses, any lifetime"),
                (r"\texttt{FCA}", r"zero-lag shape, then crossing eq.", r"$\bm v$, $\ell_1$, $\ell_2$, $\alpha$",
                 r"Gaussian pulses, any tilt and lifetime")]
        for k, row in enumerate(rows):
            y = y0 - dy * (k + 1)
            table.add(*[Tex(cell, font_size=25, color=OI_GREEN if k == 1 and c == 0 else INK
                            ).move_to([x, y, 0], aligned_edge=LEFT) for c, (cell, x) in enumerate(zip(row, xs))])
        limits = VGroup(Tex(r"Limitations", font_size=30, color=OI_RED),
                        bullets(r"the zero-lag shape assumes Gaussian pulses; the crossing equation does not",
                                r"white noise only: coherent pickup and cross-talk raise $R(\bm\Delta,0)$",
                                r"a spread of velocities broadens the shape with the lag, "
                                r"$W(\tau)=W+\Sigma_V\tau^2$",
                                r"the extrapolated zero lag needs $\mathrm{d}t$ well below the correlation time",
                                size=25, buff=0.18)
                        ).arrange(DOWN, aligned_edge=LEFT, buff=0.25).next_to(table, DOWN, buff=0.6
                                                                             ).to_edge(LEFT, buff=0.9)
        self.play(FadeIn(head), FadeIn(table))
        self.next_slide()
        self.play(FadeIn(limits, lag_ratio=0.3))
