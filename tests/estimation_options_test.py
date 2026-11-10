"""Backward compatibility of the options and the new fields: no realizations."""

import inspect

import numpy as np

import velocity_estimation as ve
import velocity_estimation.full_correlation_analysis as fca
import velocity_estimation.time_delay_estimation as tde
import velocity_estimation.two_dim_velocity_estimates as td


def test_default_options():
    eo = td.EstimationOptions()
    assert eo.method == tde.TDEMethod.CC
    assert eo.use_3point_method is True
    assert eo.cache is True
    assert eo.velocity_method == td.VelocityMethod.TDE
    assert isinstance(eo.neighbour_options, td.NeighbourOptions)
    assert isinstance(eo.cc_options, tde.CCOptions)
    assert isinstance(eo.ca_options, tde.CAOptions)
    assert isinstance(eo.ccf_options, tde.CCFitOptions)
    assert isinstance(eo.fca_options, fca.FCAOptions)
    assert eo.get_time_delay_options() is eo.cc_options


def test_fca_options_defaults():
    options = fca.FCAOptions()
    assert options.window is None
    assert options.remove_noise_spike is True
    assert options.minimum_correlation == 0.1


def test_options_set_after_construction():
    eo = td.EstimationOptions()
    eo.use_3point_method = False
    assert eo.use_3point_method is False
    eo.method = tde.TDEMethod.CA
    assert eo.get_time_delay_options() is eo.ca_options
    eo.velocity_method = td.VelocityMethod.FCA
    assert eo.velocity_method == td.VelocityMethod.FCA


def test_given_options_are_kept():
    cc_options = tde.CCOptions(cc_window=3)
    fca_options = fca.FCAOptions(window=2)
    eo = td.EstimationOptions(
        cc_options=cc_options,
        fca_options=fca_options,
        velocity_method=td.VelocityMethod.CROSSING,
    )
    assert eo.cc_options is cc_options
    assert eo.fca_options is fca_options
    assert eo.velocity_method == td.VelocityMethod.CROSSING


def test_option_instances_are_not_shared():
    a, b = td.EstimationOptions(), td.EstimationOptions()
    a.neighbour_options.max_separation = 3
    a.cc_options.cc_window = 5
    a.ca_options.distance = 1
    a.ccf_options.fit_window = 7
    a.fca_options.window = 2
    assert b.neighbour_options.max_separation == 1
    assert b.cc_options.cc_window is None
    assert b.ca_options.distance == 0
    assert b.ccf_options.fit_window == 100
    assert b.fca_options.window is None


def test_str_includes_the_new_members():
    text = str(td.EstimationOptions(velocity_method=td.VelocityMethod.FCA))
    assert "VelocityMethod.FCA" in text
    assert "Minimum Correlation: 0.1" in text


def test_pixel_data_defaults():
    pixel = td.PixelData()
    for field in ("vx", "vy", "l1", "l2", "tilt"):
        assert np.isnan(getattr(pixel, field))
    assert pixel.confidence == 0
    assert pixel.events == 0
    assert pixel.is_dead is False


def test_movie_data_getters():
    for name in ("get_vx", "get_vy", "get_confidences", "get_l1", "get_l2", "get_tilt"):
        assert callable(getattr(td.MovieData, name))


def test_signatures_keep_the_positional_arguments():
    parameters = list(inspect.signature(td.estimate_velocities_for_pixel).parameters)
    assert parameters == [
        "x",
        "y",
        "ds",
        "estimation_options",
        "tde_delegator",
        "correlation_delegator",
    ]
    parameters = inspect.signature(td.estimate_velocity_field).parameters
    assert parameters["estimation_options"].default is None


def test_package_exports():
    for name in (
        "VelocityMethod",
        "FCAOptions",
        "CorrelationDelegator",
        "spike_free_variance",
        "crossing_time",
        "shape_from_zero_lag",
        "velocity_from_crossings",
        "shape_parameters",
        "EstimationOptions",
        "estimate_velocity_field",
    ):
        assert hasattr(ve, name)
