# harmonic_imputation
# Missing data imputation by harmonic decomposition

This toolbox includes the functions to implement the missing data imputation based on harmonic decomposition of non-stationary signals. The proposed algorithm is comprised of three main steps:

1 - Initial Data Imputation: An initial missing data imputation is performed using some of the implemented imputation methods included in the toolbox. Alternatively, new methods may be include in the procedure by following the input output argument structure proposed. 

2 - Harmonic Decomposition: Harmonic amplitude and phase functions are obtained from the imputed signal. In this step, various steps are performed including STFT computing, ridge extraction, optimal harmonic number estimation and finally harmonic reconstruction is performed. 

3 - Interpolation at the Harmonic Level: Finally, the harmonic amplitude and phase functions are clipped inside the missing data interval and interpolation is performed. Two alternatives for missing data imputation are included in this implementation: cubic spline and shape-preserving cubic interpolation. 

The main functions included in this module are briefly described below. Refer to the functions documentation for details about input and output argument usage.

<code>missing_ints</code>: Automatically detect the location and length of the missing data intervals.

<code>impute</code>: Perform initial imputation using one of the implemented method or pass-by-reference your own imputation method. The passed method has to respect the input-output scheme of the implemented methods.

<code>harm_decomp</code>: Run the harmonic decomposition algorithm on the imputed signal.

<code>harm_int</code>: Interpolate at the harmonic level using one of the two possible interpolation schemes: 'spline' or 'pchip'.

<code>impute_harm_int</code>: Global function that runs all steps of the algorithm sequentially. Flags can be set to only run the initial imputation result.

## Initial Data Imputation

The first stage of the proposed algorithm includes the detection of the missing data intervals and the initial imputation performed used one of the previously implemented methods.

The missing data intervals are located using the function <code>missing_ints</code>. The input parameters of this function are the signal under study and the type of missing data to be detected. Additionally, the missing data detection can be fine-tuned by incorporating additional parameters. See the function documentation for more details.

The initial imputation is performed using the function <code>impute</code>, by passing the missing data signal and missing data intervals as parameters. A series of missing data imputation procedures are implemented:

* Takens' Lag Map (TLM).
* Dynamical System Forecasting by Least-Square Estimation (LSE).
* Dynamical Mode Decomposition (DMD).
* Extended Dynamical Mode Decomposition (EDMD).
* Gaussian Process Regression (GPR).
* ARIMA Forecasting with Forward Step (ARIMAF).
* ARIMA Forecasting with Backward Step (ARIMAB).
* Trigonometric Box-Cox, ARMA and Seasonal Modelling (TBATS).
* Sparse Time-Frequency Non-linear Matching Pursuit (STF).
* Locally Stationary Wavelet Process (LSW).

In addition, new imputation methods can be passed as inputs of the <code>impute</code> function. This methods have to follow the input/output structure of the already implemented methods. See the function documentation for more details.

## Trend Separation and Harmonic Decomposition

The initially imputed signal is decomposed using the <code>harm_decomp</code> function by computing its spectrogram and reconstructing each harmonic AM-FM component of the signal by way of spectrogram ridge detection and band-limited reconstruction. Likewise, the trend component of the signal is extracted by time-varying low-pass filtering of the imputed signal.
Parameters for the harmonic decomposition include STFT window length, reconstruction band-width and number of harmonic components K.

## Interpolation at the Harmonic Level

The harmonic amplitude and phases of the imputed signal are obtained from the extracted modes and an enhanced imputation is obtained by direct interpolation of the amplitudes and phases using the <code>harm_int</code> function. Different interpolation schemes can be applied to this task, including linear, cubic spline and cubic hermite function (pchip) interpolation. 


#IMPORTANT! Running R scripts using matlab

Some impute methods (namely, TBATS and LSW) are implemented in R. In order to run these methods, R version 2.14.0 or above has to be installed in the 'usr/bin' repository (Ubuntu) and the R.matlab ([Link](guides/content/editing-an-existing-page](https://cran.r-project.org/web/packages/R.matlab/index.html)) repository has to be installed. For more information about running R scripts in Matlab, plase refer to (https://nonlinear.wtu.edu.cn/info/1085/1055.htm).
