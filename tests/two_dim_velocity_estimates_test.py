import numpy as np
import velocity_estimation.time_delay_estimation as tde
import velocity_estimation.two_dim_velocity_estimates as td
import velocity_estimation.utils as u
import test_utils as tu
import xarray as xr


def get_estimation_options():
    return td.EstimationOptions(
        method=tde.TDEMethod.CC,
        use_3point_method=True,
    )


def test_rad_and_pol():
    v, w = 1, 1
    ds = tu.make_2d_realization(v, w, np.array([5, 6, 7]), np.array([5, 6, 7]))
    estimation_options = get_estimation_options()
    pd = td.estimate_velocities_for_pixel(
        1, 1, u.SyntheticBlobImagingDataInterface(ds), estimation_options
    )
    (
        v_est,
        w_est,
    ) = (
        pd.vx,
        pd.vy,
    )
    error = np.max([abs(v_est - v), abs(w_est - w)])
    assert error < 0.1, "Numerical error too big"


def test_full():
    v, w = 1, 1
    ds = tu.make_2d_realization(v, w, np.array([5, 6, 7]), np.array([5, 6, 7, 8]))
    eo = get_estimation_options()
    movie_data = td.estimate_velocity_field(u.SyntheticBlobImagingDataInterface(ds), eo)
    vx = movie_data.get_vx()
    assert np.max(np.abs(vx - np.ones(shape=(4, 3)))) < 0.1, "Numerical error too big"


def test_rad_and_neg_pol():
    v, w = 1, -1
    ds = tu.make_2d_realization(v, w, np.array([5, 6, 7]), np.array([5, 6, 7]))
    pd = td.estimate_velocities_for_pixel(
        1, 1, u.SyntheticBlobImagingDataInterface(ds), get_estimation_options()
    )
    (
        v_est,
        w_est,
    ) = (
        pd.vx,
        pd.vy,
    )
    error = np.max([abs(v_est - v), abs(w_est - w)])
    assert error < 0.1, "Numerical error too big"


def test_rad_and_2pol():
    v, w = 1, 2
    ds = tu.make_2d_realization(v, w, np.array([5, 6, 7]), np.array([5, 6, 7]))
    pd = td.estimate_velocities_for_pixel(
        1, 1, u.SyntheticBlobImagingDataInterface(ds), get_estimation_options()
    )
    (
        v_est,
        w_est,
    ) = (
        pd.vx,
        pd.vy,
    )
    error = np.max([abs(v_est - v), abs(w_est - w)])
    assert error < 0.1, "Numerical error too big"


def test_cond_av():
    v, w = 1, -1
    ds = tu.make_2d_realization(v, w, np.array([5, 6, 7]), np.array([5, 6, 7]))
    estimation_options = get_estimation_options()
    cond_av_eo = tde.CAOptions(length_of_return=5, distance=5, min_threshold=0.6)
    estimation_options.method = tde.TDEMethod.CA
    estimation_options.ca_options = cond_av_eo
    pd = td.estimate_velocities_for_pixel(
        1, 1, u.SyntheticBlobImagingDataInterface(ds), estimation_options
    )
    (
        v_est,
        w_est,
    ) = (
        pd.vx,
        pd.vy,
    )
    error = np.max([abs(v_est - v), abs(w_est - w)])
    assert error < 0.15, "Numerical error too big"


# Dead pixels have already been preprocessed and have an array of nans at their site
class MockXDS(xr.Dataset):
    def __init__(self, ds):
        super().__init__(ds)
        self.ds = ds

    def isel(
        self,
        indexers=None,
        drop: bool = False,
        missing_dims: str = "raise",
        **indexers_kwargs
    ):
        if indexers_kwargs["x"] == 0:
            dummy = np.array([np.nan, np.nan, np.nan])
            return xr.Dataset({"n": (["t"], dummy)})
        return self.ds.isel(indexers, drop, missing_dims, **indexers_kwargs)


def test_ignore_dead_pixels():
    v, w = 1, 1
    ds = tu.make_2d_realization(v, w, np.array([5, 6, 7]), np.array([5, 6, 7]))
    mock_ds = MockXDS(ds)
    pd = td.estimate_velocities_for_pixel(
        1, 1, u.SyntheticBlobImagingDataInterface(mock_ds), get_estimation_options()
    )
    (
        v_est,
        w_est,
    ) = (
        pd.vx,
        pd.vy,
    )
    error = np.max([abs(v_est - v), abs(w_est - w)])
    assert error < 0.1, "Numerical error too big"


def test_neighbours():
    v, w = 1.05, 1.05
    ds = tu.make_2d_realization(
        v, w, np.array([5, 5.1, 5.2]), np.array([5, 5.1, 5.2]), dt=0.1
    )
    estimation_options = get_estimation_options()
    estimation_options.neighbour_options.ccf_min_lag = 1
    estimation_options.neighbour_options.max_separation = 3
    pd = td.estimate_velocities_for_pixel(
        0, 0, u.SyntheticBlobImagingDataInterface(ds), estimation_options
    )
    (
        v_est,
        w_est,
    ) = (
        pd.vx,
        pd.vy,
    )
    error = np.max([abs(v_est - v), abs(w_est - w)])
    assert error < 0.1, "Numerical error too big"


def test_cross_corr_fit():
    v, w = 1, 1
    ds = tu.make_2d_realization(v, w, np.array([5, 6, 7]), np.array([5, 6, 7]))
    eo = get_estimation_options()
    eo.method = tde.TDEMethod.CCFit
    pd = td.estimate_velocities_for_pixel(
        1, 1, u.SyntheticBlobImagingDataInterface(ds), eo
    )
    (
        v_est,
        w_est,
    ) = (
        pd.vx,
        pd.vy,
    )
    error = np.max([abs(v_est - v), abs(w_est - w)])
    assert error < 0.1, "Numerical error too big"


def test_non_orthogonal_points():
    """Pixels on a sheared grid, (0, 5) and (1, 5) in one row, (1.5, 6) and (2.5, 6)
    in the next: picked from one realization on the fine uniform grid
    x = 0, 0.5, ..., 2.5, y = 5, 6, and handed to the pipeline in the imaging
    layout frames(y, x, time) with 2-D coordinates R, Z."""
    vx, vy = 1, 1
    x_fine = np.arange(0, 3, 0.5)
    ds = tu.make_2d_realization(vx, vy, x_fine, np.array([5.0, 6.0]))

    # Columns of x_fine holding the two pixels of each row: x = 0, 1 at y = 5 and
    # x = 1.5, 2.5 at y = 6.
    columns = np.array([[0, 2], [3, 5]])
    n = ds["n"].values
    frames = np.stack([n[row, columns[row], :] for row in range(2)])
    R = x_fine[columns]
    Z = np.repeat(ds.y.values[:, np.newaxis], 2, axis=1)

    ds_new = xr.Dataset(
        {"frames": (["y", "x", "time"], frames)},
        coords={
            "R": (["y", "x"], R),
            "Z": (["y", "x"], Z),
            "time": (["time"], ds.t.values),
        },
    )

    estimation_options = get_estimation_options()
    estimation_options.cc_options.cc_window = 1000

    movie_data = td.estimate_velocity_field(
        u.CModImagingDataInterface(ds_new), estimation_options
    )
    vx = movie_data.get_vx()
    assert np.max(np.abs(vx - np.ones(shape=(2, 2)))) < 0.1, "Numerical error too big"
