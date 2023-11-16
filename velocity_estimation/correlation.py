import numpy as np
import scipy.signal as ssi


def corr_fun(x, y, dt, norm=True, biased=True, method="auto"):
    """Estimates the correlation function between X and Y using ssi.correlate.
    For now, we require both signals to be of equal length.

    Input:
        X: First array to be correlated. ............................. (Nx1) np.array
        Y: Second array to be correlated. ............................ (Nx1) np.array
        dt: Time step of the time series. ............................ float
        norm: Normalizes the correlation function to a maxima of 1 ... bool
        biased: Trigger estimator biasing. ........................... bool
        method: 'direct', 'fft' or 'auto'. Passed to ssi.correlate ... string

    For biased=True, the result is divided by X.size.
    For biased=False, the estimator is unbiased and returns the result
    divided by X.size-|k|, where k is the lag.
    The unbiased estimator diverges for large lags, and
    for small lags and large X.size, the difference is trivial.
    """

    assert x.size == y.size
    length = x.size

    if norm:
        x = (x - x.mean()) / x.std()
        y = (y - y.mean()) / y.std()

    ccf = ssi.correlate(x, y, mode="full", method=method)

    k = np.arange(-(length - 1), length)
    times = k * dt

    ccf /= length if biased else length - np.abs(k)

    return times, ccf
