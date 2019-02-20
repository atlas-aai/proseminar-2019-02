# Psychometric Considerations for Learning Maps-Based Assessments

Example map validation analyses using simulated data. For additional resources, see []().

`01-sim-data.R` includes functions to simulate MIRT and LCDM data when using a simple Q-matrix desing.

`02-estimate-models.R` includes code for estimating MIRT, LCDM, and LCA models using [**rstan**](http://mc-stan.org/rstan/). The code for the *Stan* models can be found in `stan-models/`.

`03-map-validation.R` includes code for running the map validation analyses including posterior predictive model checks, attriubte mastery, and attribute difficulty.

`04-extra-viz.R` includes code for generating other graphics in the presentation.

All generated figures can be found in the `figures/` directory.
