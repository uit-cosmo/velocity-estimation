.. _velocity_estimation-estimation_options:

Estimation Options
===============

The class ``EstimationOptions`` allows to set the parameters for the estimation, these are:

.. list-table:: EstimationOptions
   :widths: 25 25 50
   :header-rows: 1

   * - Member
     - Type (default)
     - Description
   * - method
     - tde.TDEMethod
     - Specifies the time delay method to be used
   * - use_3point_method
     - bool(True)
     - If False, use 2-point method for velocity estimation
   * - cache
     - bool(True)
     - If True, TDE results are cached
   * - neighbour_options
     - NeighbourOptions
     - Options for neighbour selection
   * - cc_options
     - tde.CCOptions
     - Used if method = TDEMethod.CC
   * - ca_options
     - tde.CAOptions
     - Used if method = TDEMethod.CA
   * - ccf_options
     - tde.CCFitOptions
     - Used if method = TDEMethod.CCFit

As you can see, many of these options are also package classes, as ``NeighbourOptions`` or ``NeighbourOptions``. These are further explained below. Note that the member ``method`` is used to specify which method is employed for time-delay estimation (TDE), these can be based on cross-correlation (CC), conditional averaging (CA) and cross-correlation fitting (CCf), once this is set only the options pertaining the relevant method will be used (``cc_options``, ``ca_options`` or ``ccf_options``).

Neighbour Options
===============

Used for neighbour selection

.. list-table:: NeighbourOptions
   :widths: 25 25 50
   :header-rows: 1

   * - Member
     - Type (default)
     - Description
   * - ccf_min_lag
     - int(-1)
     - Checks that maximal correlation time at least ccf_min_lag time the sampling time
   * - min_separation
     - int(1)
     - Minimum allowed separation between pixels
   * - max_separation
     - int(1)
     - Maximum allowed separation

CC Options
===============

Used if the time-delay estimation is based on cross-correlations (``method=tde.TDEMethod.CC``)


.. list-table:: CCOptions
   :widths: 25 25 50
   :header-rows: 1

   * - Member
     - Type (default)
     - Description
   * - cc_window
     - float(none)
     - at Time lag window for cross-correlation; defaults to 100 samples if None
   * - minimum_cc_value
     - float(0.5)
     - Minimum value for cross-correlation maximum; no estimate if less
   * - running_mean
     - bool(True)
     - If True, applies running mean to cross-correlation
   * - running_mean_window_max
     - int(7)
     - Max length for running mean; warning if exceeded
   * - interpolate
     - bool(False)
     - If True, maximizing time lags found by interpolation


CA Options
===============

Used if the time-delay estimation is based on conditional averaging (``method=tde.TDEMethod.CA``)


.. list-table:: CAOptions
   :widths: 25 25 50
   :header-rows: 1

   * - Member
     - Type (default)
     - Description
   * - min_threshold
     - float(2.5)
     - min threshold for conditional averaged events
   * - max_threshold
     - float(np.inf)
     - max threshold for conditional averaged events
   * - delta
     - float(None)
     - If window = True, delta is the minimal distance between two peaks.
   * - window
     - bool(False)
     - If True, delta also gives the minimal distance between peaks.
   * - interpolate
     - bool(False)
     - If True the maximizing time lags are found by interpolation.
   * - verbose
     - bool(False)
     - If True prints event info


CCF Options
===============

Used if the time-delay estimation is based on cross-correlation fitting (``method=tde.TDEMethod.CCF``)


.. list-table:: CAOptions
   :widths: 25 25 50
   :header-rows: 1

   * - Member
     - Type (default)
     - Description
   * - fit_window
     - float(100)
     - The window employed for the fit will be centered at 0 with length 2 * fit_window + 1.
   * - initial_guess
     - np.array([1, 0, 1])
     - Initial guess for the optimization procedure.
   * - interpolate
     - bool(False)
     - If True the maximizing time lags are found by interpolation.
