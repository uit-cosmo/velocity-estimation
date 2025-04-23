def test_example():
    # BEGIN EXAMPLE velocity_estimation
    import numpy as np
    import velocity_estimation.time_delay_estimation as tde
    import velocity_estimation.two_dim_velocity_estimates as ve
    import velocity_estimation.utils as u
    import test_utils as tu
    import xarray as xr

    v, w = 1, 1

    # create synthetic dataset with velocities v, w evaluated on a grid of 3x3 points.
    synthetic_dataset = tu.make_2d_realization(
        v, w, np.array([5, 6, 7]), np.array([5, 6, 7])
    )

    # create a wrapper of the data, since this is synthetic data we use SyntheticBlobImagingDataInterface
    dataset_wrapper = u.SyntheticBlobImagingDataInterface(synthetic_dataset)

    # our estimates will use a cross-correlation method to estimate time delays, and it will employ the
    # three-point method.
    estimation_options = ve.EstimationOptions(
        method=tde.TDEMethod.CC,
        use_3point_method=True,
    )

    # estimate the velocity field
    velocity_field = ve.estimate_velocity_field(dataset_wrapper, estimation_options)

    # retrieve the estimated velocity components. For each component, the result is a 3x3 matrix, with one estimated
    # velocity component per pixel
    estimated_vs = velocity_field.get_vx()
    estimated_ws = velocity_field.get_vy()

    # check that the estimated velocities are close to the input velocities
    error = np.max([abs(estimated_vs.mean() - v), abs(estimated_ws.mean() - w)])
    assert error < 0.1, "Numerical error too big"
    # END EXAMPLE velocity_estimation
