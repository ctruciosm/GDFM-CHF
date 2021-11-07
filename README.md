# GDFM-CHF
These codes reproduce the empirical example in Trucíos et al. (2021)

The main file is `APP_GDFM_CHF_DCC_NL_bic.m`, where our proposal is implemented. The output is an array with the full conditional covariances matrices. Those matrices 
are the input data in `Performance_OoS.R` file (which performs the minimum variance portfolio optimization).

Alternative procedures also used in the empirical application (most of them) are in the paste `alternative_procedures`. The implementations in the manuscript started several years ago and nowadays some methods implemented in R are broken because some functions/packages are deprecated or changed. I will try to re-implement these codes but I have no time at the moment. Please, be patient.

# References

- Trucíos, C., Mazzeu, J. H., Hallin, M., Hotta, L. K., Valls Pereira, P. L., & Zevallos, M. (2021). Forecasting conditional covariance matrices in high-dimensional time series: a general dynamic factor approach. Journal of Business & Economic Statistics, (Accepted), 1-35.
