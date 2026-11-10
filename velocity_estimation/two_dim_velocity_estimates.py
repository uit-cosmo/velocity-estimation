import warnings

import velocity_estimation.time_delay_estimation as tde
import velocity_estimation.full_correlation_analysis as fca
from velocity_estimation import utils
import numpy as np
from dataclasses import dataclass
from enum import Enum
import velocity_estimation.correlation as cf


@dataclass
class NeighbourOptions:
    def __init__(
        self, ccf_min_lag: int = -1, min_separation: int = 1, max_separation: int = 1
    ):
        """Neighbour selection algorithm: For each reference pixel P0, four
        combinations of neighbouring pixels are selected to estimate the
        velocity: (up, right), (up, left), (down, right) and (down, left). For
        each combinations, the two neighbouring pixels plus the reference pixel
        are used to estimate velocities. The resulting velocity is estimated as
        the mean of all resulting velocities. If a given combination does not
        lead to a velocity estimate (for example, the cross-correlation has no
        maximum) then that combination will not contribute to the final mean.
        If the reference pixel is on the boundary, only two combinations are
        available. If the reference pixel is on a corner only one combination
        will be available.

        In the default case, each neighbour pixel (up, down, right, left) is selected as the
        nearest neighbour in that direction. This class contains options to further control
        neighbour selection.

        - ccf_min_lag: Integer, checks that the maximal correlation between adjacent
        pixels occurs at a time larger or equal than neighbors_ccf_min_lag multiples of the discretization
        time. If that's not the case, the next neighbor will be used, and so on until a
        neighbor pixel is found compliant to this condition. If set to 0, it only checks that
        the pixel is not dead.
        - min_separation: Integer, minimum allowed separation between pixels.
        - max_separation: Integer, maximum allowed separation between pixels. If some
        condition is required (such as ccf_min_lag) and not fulfilled for pixels closer
        or at than max_separation, then no neighbors will be used and the subset of pixels
        under process will yield no estimate. The condition applies on a closed interval,
        meaning that pixels separated exactly max_separation are allowed, thus setting
        min_separation = max_separation guarantees that the neighbour pixels will be separated
        by that distance.
        """
        self.ccf_min_lag = ccf_min_lag
        self.min_separation = min_separation
        self.max_separation = max_separation

    def __str__(self):
        """Return a string representation of the NeighbourOptions object."""
        return (
            f"CCF Min Lag: {self.ccf_min_lag}, "
            f"Max Separation: {self.max_separation}, "
            f"Min Separation: {self.min_separation}"
        )


class VelocityMethod(Enum):
    """Velocity estimation methods.

    TDE = time delays, 2- or 3-point according to use_3point_method,
    CROSSING = crossing times with W^-1 = I, the velocity component along each separation,
    FCA = crossing times corrected by the measured pulse shape (full correlation analysis)
    """

    TDE = 1
    CROSSING = 2
    FCA = 3


@dataclass
class EstimationOptions:
    def __init__(
        self,
        method: tde.TDEMethod = tde.TDEMethod.CC,
        use_3point_method: bool = True,
        cache: bool = True,
        neighbour_options: NeighbourOptions = None,
        cc_options: tde.CCOptions = None,
        ca_options: tde.CAOptions = None,
        ccf_options: tde.CCFitOptions = None,
        velocity_method: VelocityMethod = VelocityMethod.TDE,
        fca_options: fca.FCAOptions = None,
    ):
        """Estimation options for velocity estimation method.

        - method: fppanalysis.time_delay_estimation.TDEMethod Specifies the time delay method to be used.
        - use_3point_method: [bool] If False, use 2 point method to estimate velocities from time delays.

        - cache: bool, if True TDE results and correlations are cached
        - neighbour_options: NeighbourOptions Options used for neighbour selection
        - cc_options: Cross correlation estimation options to be used if method = TDEMethod.CC
        - ca_options: Conditional average estimation options to be used if method = TDEMethod.CA
        - ccf_options: Time delay estimation options to be used if method = TDEMethod.CCFit
        - velocity_method: VelocityMethod, time delays (TDE, the default), crossing times (CROSSING) or
        full correlation analysis (FCA). method, use_3point_method and the time delay options are only
        used by TDE.
        - fca_options: FCAOptions Options of the CROSSING and FCA methods

        Option arguments left as None are built with their defaults, one instance per EstimationOptions.
        """
        self.method = method
        self.use_3point_method = use_3point_method
        self.cache = cache
        self.neighbour_options = (
            neighbour_options if neighbour_options is not None else NeighbourOptions()
        )
        self.cc_options = cc_options if cc_options is not None else tde.CCOptions()
        self.ca_options = ca_options if ca_options is not None else tde.CAOptions()
        self.ccf_options = (
            ccf_options if ccf_options is not None else tde.CCFitOptions()
        )
        self.velocity_method = velocity_method
        self.fca_options = fca_options if fca_options is not None else fca.FCAOptions()

    def get_time_delay_options(self):
        match self.method:
            case tde.TDEMethod.CC:
                return self.cc_options
            case tde.TDEMethod.CA:
                return self.ca_options
            case tde.TDEMethod.CCFit:
                return self.ccf_options

    def __str__(self):
        """Return a string representation of the EstimationOptions object."""
        return (
            f"Velocity Method: {self.velocity_method}, "
            f"Method: {self.method}, "
            f"Use 3-Point Method: {self.use_3point_method}, "
            f"Cache: {self.cache}, "
            f"Neighbor Options: {self.neighbour_options}, "
            f"CC Options: {str(self.cc_options)}, "
            f"CA Options: {str(self.ca_options)}, "
            f"CCF Options: {str(self.ccf_options)}, "
            f"FCA Options: {str(self.fca_options)}"
        )


@dataclass
class PixelData:
    """Data class containing estimated data from a single pixel.

    vx: Radial velocities
    vy: Poloidal velocities
    confidences:
        if method='cross_corr', CROSSING or FCA:
            Maximum value of the cross-correlations at each pixel.
        if method='cond_av':
            Conditional variance value at maximum cross conditional average for each pixel.
    R: Radial positions
    Z: Poloidal positions
    is_dead: True if pixel is dead
    l1, l2: Major and minor widths of the pulse autocorrelation, FCA only (NaN otherwise)
    tilt: Angle of its long axis from the radial coordinate, radians in (-pi/2, pi/2], FCA only
    """

    r_pos: float = 0
    z_pos: float = 0
    vx: float = np.nan
    vy: float = np.nan
    confidence: float = 0
    events: int = 0
    is_dead: bool = False
    l1: float = np.nan
    l2: float = np.nan
    tilt: float = np.nan


class MovieData:
    """Class containing estimated data for all pixels in a set. Return object
    from estimate_velocity_field function, the indexing format of MovieData
    matches that of the dataset input of estimate_velocity_field.

    Use getters to retrieve:
        vx: Radial velocities
        vy: Poloidal velocities
        confidences:
            if method='cross_corr':
                Maximum value of the cross-correlations at each pixel.
            if method='cond_av':
                Conditional variance value at maximum cross conditional average for each pixel.
        R: Radial positions
        Z: Poloidal positions
        is_dead: True if pixel is dead
        l1, l2, tilt: Widths and tilt of the pulse autocorrelation, FCA only (NaN otherwise)

    Dead pixels have empty PixelData (null vx and vy).
    """

    def __init__(
        self, ds: utils.ImagingDataInterface, estimation_options: EstimationOptions
    ):
        len_r, len_z = ds.get_shape()
        range_r, range_z = range(0, len_r), range(0, len_z)
        self.ds = ds
        self.estimation_options = estimation_options
        self.tde_delegator = tde.TDEDelegator(
            estimation_options.method,
            estimation_options.get_time_delay_options(),
            estimation_options.cache,
        )
        self.correlation_delegator = fca.CorrelationDelegator(
            estimation_options.fca_options, estimation_options.cache
        )
        self.pixels = [[PixelData() for _ in range_r] for _ in range_z]

        for i in range_z:
            for j in range_r:
                self.pixels[i][j] = self._set_pixel((j, i))

    def _set_pixel(self, items):
        i, j = items[0], items[1]
        try:
            return estimate_velocities_for_pixel(
                i,
                j,
                self.ds,
                self.estimation_options,
                self.tde_delegator,
                self.correlation_delegator,
            )
        except:
            print(
                "Issues estimating velocity for pixel",
                i,
                j,
                "Run estimate_velocities_for_pixel(i, j, ds, eo) to get a detailed error stacktrace",
            )
        return PixelData()

    def _get_field(self, field_name):
        return np.array(
            [[getattr(p, field_name) for p in pixel_row] for pixel_row in self.pixels]
        )

    def get_vx(self):
        return self._get_field("vx")

    def get_vy(self):
        return self._get_field("vy")

    def get_R(self):
        return self._get_field("r_pos")

    def get_Z(self):
        return self._get_field("z_pos")

    def get_events(self):
        return self._get_field("events")

    def get_confidences(self):
        return self._get_field("confidence")

    def get_is_dead(self):
        return self._get_field("is_dead")

    def get_l1(self):
        return self._get_field("l1")

    def get_l2(self):
        return self._get_field("l2")

    def get_tilt(self):
        return self._get_field("tilt")


def get_2d_velocities_from_time_delays(
    delta_t1, delta_t2, delta_x1, delta_y1, delta_x2, delta_y2
):
    """
    Estimates radial and poloidal velocities given the input parameters:
    Input:
         delta_t1 Estimation of the time delay between radially separated points.
         delta_t2 Estimation of the time delay between poloidally separated points.
         delta_x1 Radial separation between radially separated points.
         delta_y1 Poloidal separation between radially separated points.
         delta_x2 Spatial separation between poloidally separated points.
         delta_y2 Poloidal separation between poloidally separated points.

    These quantities should be obtained from three pixel points: a reference pixel point,
    a pixel point separated radially, and a pixel point separated poloidally.
    Returns:
         vx Radial velocity
         vy Poloidal velocity
    """

    numerator_v_term = delta_t2 * delta_y1 - delta_t1 * delta_y2
    numerator_w_term = delta_t2 * delta_x1 - delta_t1 * delta_x2
    common_denominator_term = numerator_w_term**2 + numerator_v_term**2
    determinant = delta_x1 * delta_y2 - delta_x2 * delta_y1

    # Ensure the common denominator term is not zero to avoid division by zero error
    if common_denominator_term == 0:
        raise ValueError(
            "The common denominator term is zero, which will lead to a division by zero error."
        )

    if determinant == 0:
        warnings.warn(
            "Determinant of separations vanishes! This will result in vanishing velocities"
        )

    v = -numerator_v_term * determinant / common_denominator_term
    w = numerator_w_term * determinant / common_denominator_term

    return v, w


def get_1d_velocities_from_time_delays(delta_tx, delta_ty, delta_x, delta_y):
    """
    Estimates radial and poloidal velocities from naive method
    given the input parameters:
    Input:
         delta_tx Estimation of the time delay between radially separated points.
         delta_ty Estimation of the time delay between poloidally separated points.
         delta_x Spatial separation between radially separated points.
         delta_y Spatial separation between poloidally separated points.

    These quantities should be obtained from two pixel points:
        radial direction: a reference pixel point and a pixel point separated radially
        poloidal direction: a reference pixel point and a pixel point separated poloidally.
    Returns:
         vx Radial velocity
         vy Poloidal velocity
    """
    vx = 0 if delta_tx == 0 else delta_x / delta_tx
    vy = 0 if delta_ty == 0 else delta_y / delta_ty

    return vx, vy


def _estimate_velocities_given_points(
    p0,
    p1,
    p2,
    ds: utils.ImagingDataInterface,
    tde_delegator: tde.TDEDelegator,
    use_2d_estimation: bool,
):
    """Estimates radial and poloidal velocity from estimated time delay either
    from cross conditional average between the pixels or cross correlation.

    This is specified in method argument.
    """
    delta_ty, cy, events_y = tde_delegator.estimate_time_delay(p2, p0, ds)
    delta_tx, cx, events_x = tde_delegator.estimate_time_delay(p1, p0, ds)

    # If for some reason the time delay cannot be estimated, we return None
    if delta_tx is None or delta_ty is None:
        return None

    confidence = min(cx, cy)
    events = min(events_x, events_y)

    r0, z0 = ds.get_position(p0[0], p0[1])
    r1, z1 = ds.get_position(p1[0], p1[1])
    r2, z2 = ds.get_position(p2[0], p2[1])

    if use_2d_estimation:
        return (
            *get_2d_velocities_from_time_delays(
                delta_tx, delta_ty, r1 - r0, z1 - z0, r2 - r0, z2 - z0
            ),
            confidence,
            events,
        )
    else:
        return (
            *get_1d_velocities_from_time_delays(delta_tx, delta_ty, r1 - r0, z2 - z0),
            confidence,
            events,
        )


def _check_ccf_constrains(
    p0, p1, ds: utils.ImagingDataInterface, neighbors_ccf_min_lag: int
):
    """Returns true if the time lag that maximizes the cross-correlation
    function measure at p0 and p1 is not zero."""
    if ds.is_pixel_dead(p1[0], p1[1]):
        return False

    signal0 = ds.get_signal(p0[0], p0[1])
    signal1 = ds.get_signal(p1[0], p1[1])

    # No need to compute the ccf if the min lag is 0
    if neighbors_ccf_min_lag == 0:
        return True

    ccf_times, ccf = cf.corr_fun(
        signal0, signal1, dt=ds.get_dt(), biased=True, norm=True
    )
    ccf = ccf[np.abs(ccf_times) < max(ccf_times) / 2]
    ccf_times = ccf_times[np.abs(ccf_times) < max(ccf_times) / 2]
    max_index = np.argmax(ccf)

    fulfills_constrain = (
        np.abs(ccf_times[max_index]) >= neighbors_ccf_min_lag * ds.get_dt()
    )
    if not fulfills_constrain:
        warnings.warn(
            "Pixel {} does not fulfill cross-correlation time lag condition with respect to pixel of {}."
            " Updating.".format(p1, p0)
        )

    return fulfills_constrain


def _find_neighbors(
    x, y, ds: utils.ImagingDataInterface, neighbour_options: NeighbourOptions
):
    start = neighbour_options.min_separation
    end = neighbour_options.max_separation

    def fulfills_conditions(p):
        return ds.is_within_boundaries(p[0], p[1]) and _check_ccf_constrains(
            (x, y), p, ds, neighbour_options.ccf_min_lag
        )

    horizontal = []
    vertical = []

    left = -start
    while np.abs(left) <= end:
        if fulfills_conditions((x + left, y)):
            horizontal.append((x + left, y))
        left -= 1

    right = start
    while np.abs(right) <= end:
        if fulfills_conditions((x + right, y)):
            horizontal.append((x + right, y))
        right += 1

    up = start
    while np.abs(up) <= end:
        if fulfills_conditions((x, y + up)):
            vertical.append((x, y + up))
        up += 1

    down = -start
    while np.abs(down) <= end:
        if fulfills_conditions((x, y + down)):
            vertical.append((x, y + down))
        down -= 1

    return horizontal, vertical


def _estimate_velocities_from_crossings(
    x,
    y,
    h_neighbors,
    v_neighbors,
    ds: utils.ImagingDataInterface,
    estimation_options: EstimationOptions,
    correlation_delegator: fca.CorrelationDelegator,
):
    """Crossing time and FCA estimates for the reference pixel (x, y): one pooled fit over its
    neighbours. The crossing time of the reference with each neighbour gives one equation for
    the velocity. For the FCA, the zero-lag correlations of the reference with each neighbour
    and of each horizontal neighbour with each vertical one give the shape W^-1 that corrects
    the crossing equations; the crossing method uses the identity instead. Both a horizontal
    and a vertical neighbour are required."""
    r_pos, z_pos = ds.get_position(x, y)
    if len(h_neighbors) == 0 or len(v_neighbors) == 0:
        return PixelData(r_pos=r_pos, z_pos=z_pos)

    def separation(p1, p0):
        r1, z1 = ds.get_position(p1[0], p1[1])
        r0, z0 = ds.get_position(p0[0], p0[1])
        return float(r1) - float(r0), float(z1) - float(z0)

    p0 = (x, y)
    neighbours = h_neighbors + v_neighbors
    lags, r_auto_0 = correlation_delegator.autocorrelation(p0, ds)
    separations, crossing_times, maxima = [], [], []
    for p1 in neighbours:
        _, r_auto_1 = correlation_delegator.autocorrelation(p1, ds)
        _, r_cross = correlation_delegator.cross_correlation(p1, p0, ds)
        tau = fca.crossing_time(lags, (r_auto_0 + r_auto_1) / 2, r_cross)
        if np.isnan(tau):
            warnings.warn(
                "Cross-correlation not below the autocorrelation at zero lag, no crossing time"
                " between pixels {} and {}".format(p1, p0)
            )
            continue
        separations.append(separation(p1, p0))
        crossing_times.append(tau)
        maxima.append(np.max(r_cross))

    if estimation_options.velocity_method == VelocityMethod.FCA:
        pairs = [(p1, p0) for p1 in neighbours]
        pairs += [(pv, ph) for ph in h_neighbors for pv in v_neighbors]
        shape_separations, r0 = [], []
        for p1, p2 in pairs:
            r = correlation_delegator.zero_lag_correlation(p1, p2, ds)
            if r >= estimation_options.fca_options.minimum_correlation:
                shape_separations.append(separation(p1, p2))
                r0.append(r)
        w_inv = fca.shape_from_zero_lag(shape_separations, r0)
        if w_inv is None:
            warnings.warn(
                f"Zero-lag correlations do not determine a positive definite shape for pixel {p0}"
            )
            return PixelData(r_pos=r_pos, z_pos=z_pos)
    else:
        w_inv = np.eye(2)

    v = fca.velocity_from_crossings(separations, crossing_times, w_inv)
    if v is None:
        warnings.warn(
            "Crossing times do not determine the velocity for pixel {}".format(p0)
        )
        return PixelData(r_pos=r_pos, z_pos=z_pos)

    pixel_data = PixelData(
        r_pos=r_pos,
        z_pos=z_pos,
        vx=float(v[0]),
        vy=float(v[1]),
        confidence=float(min(maxima)),
        events=0,
    )
    if estimation_options.velocity_method == VelocityMethod.FCA:
        pixel_data.l1, pixel_data.l2, pixel_data.tilt = fca.shape_parameters(w_inv)
    return pixel_data


def estimate_velocities_for_pixel(
    x,
    y,
    ds: utils.ImagingDataInterface,
    estimation_options: EstimationOptions = None,
    tde_delegator: tde.TDEDelegator = None,
    correlation_delegator: fca.CorrelationDelegator = None,
):
    """Estimates radial and poloidal velocity for a pixel with indexes x,y
    using all four possible combinations of nearest neighbour pixels (x-1, y),
    (x, y+1), (x+1, y) and (x, y-1). Dead-pixels (stored as np.nan arrays) are
    ignored. Pixels outside the coordinate domain are ignored. Time delay
    estimation is performed by maximizing either the cross- correlation
    function or cross conditional average function, which is specified in input
    argument 'method'.

    If time delay estimation is performed by maximizing the cross correlation function,
    the confidence of the estimation is a value in the interval (0, 1) given by the
    mean of the confidences for each combination, which is given by the minimum
    of the maximums of the two cross-correlations involved.

    If time delay estimation is performed by maximizing the cross conditional average function,
    the confidence of the estimation is a value in the interval (0, 1) given by the
    cross conditional variance for each event. OBS: We return 1-CV for cross conditional variance.

    With velocity_method CROSSING or FCA, the velocity comes instead from one pooled fit of the
    crossing times of the reference with its neighbours, corrected by the pulse shape measured
    from the zero-lag correlations in the FCA case (see full_correlation_analysis). The confidence
    is the smallest maximum of the cross-correlations involved, and the FCA fills in l1, l2 and tilt.

    Input:
        x: pixel index x
        y: pixel index y
        ds: xarray Dataset
        estimation_options: EstimationOptions class including all estimation parameters, if not set
        the method will be based on cross-correlation function.
        tde_delegator, correlation_delegator: caches shared between pixels, built here if not given.


    Returns:
        PixelData: Object containing radial and poloidal velocities and method-specific data.
    """
    if estimation_options is None:
        estimation_options = EstimationOptions()

    r_pos, z_pos = ds.get_position(x, y)

    # If the reference pixel is dead, return empty data right away
    if ds.is_pixel_dead(x, y):
        return PixelData(r_pos=r_pos, z_pos=z_pos, is_dead=True)

    h_neighbors, v_neighbors = _find_neighbors(
        x, y, ds, estimation_options.neighbour_options
    )

    if estimation_options.velocity_method in (
        VelocityMethod.CROSSING,
        VelocityMethod.FCA,
    ):
        if correlation_delegator is None:
            correlation_delegator = fca.CorrelationDelegator(
                estimation_options.fca_options, estimation_options.cache
            )
        return _estimate_velocities_from_crossings(
            x,
            y,
            h_neighbors,
            v_neighbors,
            ds,
            estimation_options,
            correlation_delegator,
        )

    if tde_delegator is None:
        tde_delegator = tde.TDEDelegator(
            estimation_options.method,
            estimation_options.get_time_delay_options(),
            estimation_options.cache,
        )

    results = [
        _estimate_velocities_given_points(
            (x, y), px, py, ds, tde_delegator, estimation_options.use_3point_method
        )
        for px in h_neighbors
        if ds.is_within_boundaries(px[0], px[1])
        for py in v_neighbors
        if ds.is_within_boundaries(py[0], py[1])
    ]

    results = [r for r in results if r is not None]
    if len(results) == 0:  # If no neighbor pixels are found we cannot estimate
        return PixelData(r_pos=r_pos, z_pos=z_pos)

    mean_vx = sum(map(lambda r: r[0], results)) / len(results)
    mean_vy = sum(map(lambda r: r[1], results)) / len(results)
    confidence = sum(map(lambda r: r[2], results)) / len(results)
    events = sum(map(lambda r: r[3], results)) / len(results)

    return PixelData(
        r_pos=r_pos,
        z_pos=z_pos,
        vx=mean_vx,
        vy=mean_vy,
        confidence=confidence,
        events=events,
    )


def estimate_velocity_field(
    ds: utils.ImagingDataInterface,
    estimation_options: EstimationOptions = None,
) -> MovieData:
    """Computes the velocity field of a given dataset ds implementing a
    utils.ImagingDataInterface class. The estimation takes into account
    poloidal flows as described in the 2D filament model. For each pixel, the
    velocities are estimated using the given pixel, and two neighbour pixels:
    the right neighbour and the down neighbour. The velocities are estimated
    from a time delay estimation performed by maximizing either the cross-
    correlation function or cross conditional average function, which is
    specified in input argument 'method', or, with velocity_method CROSSING or
    FCA, from the crossing times of the correlations (see
    estimate_velocities_for_pixel).

    If time delay estimation is performed by maximizing the cross correlation function,
    the confidence of the estimation is a value in the interval (0, 1) given by the
    mean of the confidences for each combination, which is given by the minimum
    of the maximums of the two cross-correlations involved.

    If time delay estimation is performed by maximizing the cross conditional average function,
    the confidence of the estimation is a value in the interval (0, 1) given by the
    cross conditional variance for each event. OBS: We return 1-CV for cross conditional variance.

    The return objects are matrices of the size of the GPI grid,
    from which the velocity field can be easily plotted via f.e matplotlib.quiver.

    Input:
        ds: utils.ImagingDataInterface data
        estimation_options: EstimationOptions class including all estimation parameters, if not set
        the method will be based on cross-correlation function.


    Returns:
        movie_data: Class containing estimation data about all pixels
    """
    if estimation_options is None:
        estimation_options = EstimationOptions()

    return MovieData(ds, estimation_options)
