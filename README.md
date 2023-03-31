# harmonic_imputation
# Missing data imputation by harmonic decomposition

This toolbox includes the functions to implement the missing data imputation based on harmonic decomposition of non-stationary signals. The proposed algorithm is comprised of three main steps:

1 - Initial Data Imputation: An initial missing data imputation is performed using some of the implemented imputation methods included in the toolbox. Alternatively, new methods may be include in the procedure by following the input output argument structure proposed. 

2 - Harmonic Decomposition: Harmonic amplitude and phase functions are obtained from the imputed signal. In this step, various steps are performed including STFT computing, ridge extraction, optimal harmonic number estimation and finally harmonic reconstruction is performed. 

3 - Interpolation at the Harmonic Level: Finally, the harmonic amplitude and phase functions are clipped inside the missing data interval and interpolation is performed. Two alternatives for missing data imputation are included in this implementation: cubic spline and shape-preserving cubic interpolation. 

The main functions included in this module are briefly described below. Refer to the functions documentation for details about input and output argument usage.

<code>missing_ints</code>: Automatically detect the location and length of the missing data intervals.

<code>impute</code> : Perform initial imputation using one of the implemented method or pass-by-reference your own imputation method. The passed method has to respect the input-output scheme of the implemented methods.

<code>harm_decomp</code> : Run the harmonic decomposition algorithm on the imputed signal.

<code>harm_int</code> : Interpolate at the harmonic level using one of the two possible interpolation schemes: 'spline' or 'pchip'.

<code>impute_harm_int</code> : Global functions that runs all steps of the algorithm sequentially. Flags can be set to only run the initial imputation result.






