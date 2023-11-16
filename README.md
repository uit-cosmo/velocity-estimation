# velocity-estimation
Two dimensional velocity estimation methods for coarse-grained imaging data. Traditional methods compute the velocity
 components in a given directions from the time delay between signals separated in such direction. This approach is inaccurate and can lead to big errors if
  the velocity of propagation is not aligned with the separation between the two measurement points. At least
  three points need to be considered, and time delays in two different directions need to be used simultaneously for
  the accurate estimation of the velocity vector. The code in this repository implements such method for imaging
  data. The underlying time delay estimation can be switched from cross-correlation analysis or cross-conditional
  averaging.

# Install

```
git clone https://github.com/uit-cosmo/velocity-estimation.git
cd velocity-estimation
pip install .
```

# Use

The main function is two_dim_velocity_estimates.estimate_velocity_field(ds, eo), its usage is described in a notebook under the guides/ folder.
