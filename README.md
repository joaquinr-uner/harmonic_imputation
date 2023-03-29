# harmonic_imputation
# Missing data imputation by harmonic decomposition

This toolbox includes the functions to implement the missing data imputation based on harmonic decomposition of non-stationary signals. The proposed algorithm is comprised of three main steps:

1 - Initial Data Imputation: An initial missing data imputation is performed using some of the implemented imputation methods included in the toolbox. Alternatively, new methods may be include in the procedure by following the input output argument structure proposed. 

2 - Harmonic Decomposition: Harmonic amplitude and phase functions are obtained from the imputed signal. In this step, various steps are performed including STFT computing, ridge extraction, optimal harmonic number estimation and finally harmonic reconstruction is performed. 

3 - Interpolation at the Harmonic Level: Finally, the harmonic amplitude and phase functions are clipped inside the missing data interval and interpolation is performed. Two alternatives for missing data imputation are included in this implementation: cubic spline and shape-preserving cubic interpolation. 

The main functions included in this module are briefly described below. Refer to the functions documentation for details about input and output argument usage.

missing_ints :

impute : 

harm_decomp : 

harm_int : 

impute_harm_int : 






