.. _velocity_estimation-estimation_options:

Estimation Options
==================

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
   * - velocity_method
     - VelocityMethod(TDE)
     - Which estimator is used: time delays (TDE), crossing times (CROSSING) or full correlation analysis (FCA); see :ref:`velocity_estimation-estimation_options-velocity_method`
   * - fca_options
     - FCAOptions
     - Used if velocity_method = VelocityMethod.CROSSING or VelocityMethod.FCA

As you can see, many of these options are also package classes, as ``NeighbourOptions`` or ``CCOptions``. These are further explained below. Note that the member ``method`` is used to specify which method is employed for time-delay estimation (TDE), these can be based on cross-correlation (CC), conditional averaging (CA) and cross-correlation fitting (CCFit), once this is set only the options pertaining the relevant method will be used (``cc_options``, ``ca_options`` or ``ccf_options``). The member ``velocity_method`` selects the estimator itself, see :ref:`velocity_estimation-estimation_options-velocity_method`. ``method``, ``use_3point_method`` and the TDE option classes apply only to ``VelocityMethod.TDE``; ``VelocityMethod.CROSSING`` and ``VelocityMethod.FCA`` ignore them and use ``fca_options`` instead.

.. _velocity_estimation-estimation_options-neighbour:

Neighbour Options
=================

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

.. _velocity_estimation-estimation_options-cc:

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


.. _velocity_estimation-estimation_options-ca:

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
   * - length_of_return
     - float(None)
     - float, length of recorded events. If `None`, estimated as ``len(signal) / len(number_of_events) * time_step``.
   * - distance
     - float(0)
     - minimum distance between peaks
   * - interpolate
     - bool(False)
     - If True the maximizing time lags are found by interpolation.


.. _velocity_estimation-estimation_options-ccf:

CCF Options
===============

Used if the time-delay estimation is based on cross-correlation fitting (``method=tde.TDEMethod.CCF``)


.. list-table:: CCFitOptions
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


.. _velocity_estimation-estimation_options-velocity_method:

Velocity Method
===============

``velocity_method`` selects the estimator. The default, ``VelocityMethod.TDE``, is the time-delay method described on the main page, so existing scripts are unchanged. The other two are described under :ref:`velocity_estimation-biases`.

.. list-table:: VelocityMethod
   :widths: 20 50 30
   :header-rows: 1

   * - Value
     - Estimate
     - Biases removed
   * - TDE
     - The velocity from the time delays between the reference pixel and a horizontal and a vertical neighbour, estimated by ``method``, one solve per triangle, averaged over the triangles; the two-point estimate if ``use_3point_method`` is False. The default.
     - Misalignment (with ``use_3point_method``).
   * - CROSSING
     - The velocity component along each separation from the crossing time :math:`\tau_\mathrm{e}`, the lag at which the cross-correlation of the pair equals the mean of the two autocorrelations: :math:`\hat{\boldsymbol{e}}\cdot\boldsymbol{v} = \Delta/2\tau_\mathrm{e}`. The components along the separations to all neighbours are combined in one least-squares fit per reference pixel.
     - Misalignment and decorrelation. Exact for round pulses of any shape and lifetime.
   * - FCA
     - The crossing times corrected by the shape of the zero-lag cross-correlation, :math:`W^{-1}` from :math:`-2\ln R(\boldsymbol{\Delta}, 0) = \boldsymbol{\Delta}^\mathsf{T} W^{-1}\boldsymbol{\Delta}`, fitted over the pairs of the reference pixel with its neighbours and of each horizontal with each vertical neighbour; the crossing equations are then solved for the velocity in one pooled least-squares fit per reference pixel. Also returns the shape (``l1``, ``l2``, ``tilt``).
     - Misalignment, decorrelation and barber pole.

``CROSSING`` and ``FCA`` share the neighbour selection with the TDE and need both a horizontal and a vertical neighbour, as the three-point method does. They read their options from ``fca_options``; ``method``, ``use_3point_method``, ``cc_options``, ``ca_options`` and ``ccf_options`` are ignored.

.. _velocity_estimation-estimation_options-fca:

FCA Options
===============

Used if ``velocity_method`` is ``VelocityMethod.CROSSING`` or ``VelocityMethod.FCA``.


.. list-table:: FCAOptions
   :widths: 25 25 50
   :header-rows: 1

   * - Member
     - Type (default)
     - Description
   * - window
     - float(None)
     - Time lags in [-window, window], in time units, over which the autocorrelations and cross-correlations are computed; 100 sampling times if None, as ``cc_window``
   * - remove_noise_spike
     - bool(True)
     - If True, each pixel's variance is replaced by its autocovariance extrapolated to zero lag, :math:`C^*(0) = [4C(\mathrm{d}t) - C(2\,\mathrm{d}t)]/3`, and all correlations are normalized with it. This removes the zero-lag spike of white noise, which otherwise lowers every :math:`R(\boldsymbol{\Delta}, 0)`, biases the shape and shifts the crossing of two pixels with unequal noise
   * - minimum_correlation
     - float(0.1)
     - Pairs whose zero-lag cross-correlation is below this value are left out of the shape fit (FCA only)
