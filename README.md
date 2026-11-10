# velocity-estimation
Two dimensional velocity estimation methods for coarse-grained imaging data. Traditional methods compute the velocity
 components in a given directions from the time delay between signals separated in such direction. This approach is inaccurate and can lead to big errors if
  the velocity of propagation is not aligned with the separation between the two measurement points. At least
  three points need to be considered, and time delays in two different directions need to be used simultaneously for
  the accurate estimation of the velocity vector. The code in this repository implements such method for imaging
  data. The underlying time delay estimation can be switched from cross-correlation analysis or cross-conditional
  averaging.

For a superposition of pulses with a finite lifetime, elongated and tilted, the time delays are biased three ways:
misalignment of the motion with the separation, which the three-point method removes; decorrelation, a finite lifetime
pulling the maximum of the cross-correlation towards shorter lags, so that both the two- and the three-point estimates
read too high; and the barber pole, tilted structures turning the three-point estimate towards their short axis.
Version 1.1 adds two estimators. The crossing time, the lag at which the cross-correlation of a pair of pixels equals
the autocorrelation, removes the first two. Full correlation analysis (FCA) corrects the crossing times by the shape of
the zero-lag cross-correlation, removes all three and returns the shape; it follows Briggs, Phillips and Shinn, Proc.
Phys. Soc. B 63, 106 (1950), adapted to the pulse model. The three-point time-delay method stays the default; with
`import velocity_estimation as ve`, FCA is selected by

```python
eo = ve.EstimationOptions(velocity_method=ve.VelocityMethod.FCA)
movie_data = ve.estimate_velocity_field(ds, eo)
```

Documentation [here](https://velocity-estimation.readthedocs.io/en/latest/).

# Install

```
git clone https://github.com/uit-cosmo/velocity-estimation.git
cd velocity-estimation
pip install .
```

# Use

The main function is two_dim_velocity_estimates.estimate_velocity_field(ds, eo), its usage is described in a notebook under the guides/ folder.
