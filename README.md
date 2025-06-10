# PhD MATLAB code for projects associated with the 2015 Dissertation titled:

[Picosecond Time-Resolved Studies of Multiple Energy Transfer in Conjugated Polymer Nanoparticles](https://core.ac.uk/download/268651601.pdf)

## Related Journal Articles:

[Journal of Physical Chemistry C, 2017, "Effect of Swelling on Multiple Energy Transfer in Conjugated Polymer Nanoparticles"](https://pubs.acs.org/doi/10.1021/acs.jpcc.7b00892)

[Journal of Physical Chemistry C, 2013, "Measurement of Exciton Transport in Conjugated Polymer Nanoparticles"](https://pubs.acs.org/doi/10.1021/jp407065h)

## MATLAB scripts:

[Exciton Diffusion Model Notes](https://github.com/lcgroff2/matlab/blob/master/etdiffnpv3/diffmodel-notes.pdf)

Exciton Diffusion Model MATLAB scripts required to run Monte Carlo Bootstrap Poisson Statistical Model

[etmulti.m](https://github.com/lcgroff2/matlab/blob/master/etdiffnpv3/etmulti.m), [etdiffnp3.m](https://github.com/lcgroff2/matlab/blob/master/etdiffnpv3/etdiffnp3.m), [etdiffnp_fun3.m](https://github.com/lcgroff2/matlab/blob/master/etdiffnpv3/etdiffnp_fun3.m)

Atomic Force Microscopy (AFM) Image Analysis Scripts:

[AFM image processing script](https://github.com/lcgroff2/matlab/blob/master/afmpick.m), [AFM data loading script](https://github.com/lcgroff2/matlab/blob/master/loadafm.m)

Time-Correlated Single-Photon Counting (TCSPC):

The input parameters for the picosecond TCSPC fitting code is given via the link below. To use, we require:
```matlab

% User must edit the following lines.

%An ASCII data file and IRF file:
DataFile='tcspc_measured_in.asc'; % don't get the data and irf files
                                  % mixed up, or strange results occur.
IRFFile='irf_measurement_in.asc'; % to get more fit results and possibly
                                  % better statistics, try using two
                                  % different IRF scans (the one before
                                  -% the sample scan and the one after)

dt=2.6333;                % experimental dwell time in ps (picoseconds per
                          % bin). This should be re-measured
                          % periodically.

SaveFile = [DataFile '.fit.mat']; % where to save fit results.
dran=1580:1850; % Data range.  Carefully edit to include 100-300 ps
                % before t0, and to not go much past 4*tau.  Try
                % adjusting dran to test fit: if the fit is robust,
                % changing dran by 50 bins or so on either end
                % should only change the time constant by a few 
                % percent or less.
Rebinlvl=0; % for re-binning data (boxcar), esp. useful for
            % increasing the SNR for samples with long lifetimes.
            % Set Rebinlvl=0 for no re-binning.  Set to 5 for x5
            % reduction in the number of data points.
DoErrorbar=0;  % not implemented yet.  For now, get your error bars
               % by fitting several different sample and IRF runs and
               % manually calculating the standard deviation for a series
               % of sample measures.
clear fitinfo
curfun=0;

% the "t0 shift" section -- sometimes the IRF needs to be shifted
% by several ps.  Improper t0 shift leads to a sharp feature in the
% residuals that looks like the derivative of the IRF
fitinfo.t0lbound=-60;
fitinfo.t0ubound=60;
fitinfo.t0curval=-20;

% % exponential section -- uncomment following block of code to fit an
% % exponential, and edit the bounds and curval (initial guess)
curfun=curfun+1;
fitinfo.fun{curfun}.function = 'exp'; % example of exponential
fitinfo.fun{curfun}.lbound=[7];
fitinfo.fun{curfun}.ubound=[100];
fitinfo.fun{curfun}.curval=[20];
% 
% % second exponential (bi-exponential) section -- 
% % uncomment following block of code add fit a second
% % exponential, and edit the bounds and curval (initial guess)
% % BE CAREFUL THAT THE BOUNDS AND GUESS DO NOT OVERLAP WITH THOSE OF
% % THE FIRST EXPONENTIAL, OTHERWISE DIVISION BY ZERO CAN RESULT!!!
curfun=curfun+1;
fitinfo.fun{curfun}.function = 'exp'; % example of exponential
fitinfo.fun{curfun}.lbound=[101];
fitinfo.fun{curfun}.ubound=[1500];
fitinfo.fun{curfun}.curval=[350];

% kww section -- uncomment following block of code to fit a
% kww, and edit the bounds and curval (initial guess)
% curfun=curfun+1;
% fitinfo.fun{curfun}.function = 'kww'; % can be 'exp', 'kww', 'offset',
%                                       %'linear', 'hside', or 'irf'
% fitinfo.fun{curfun}.lbound=[3 .1]; % first value is the tau
%                                      %parameter, second is beta
% fitinfo.fun{curfun}.ubound=[300 1.0]; % first value is the tau
%                                        %parameter, second is beta
% fitinfo.fun{curfun}.curval=[3 .2]; % first value is the tau
%                                      %parameter, second is beta
```

The fit routine works by truncating the data to within the data region (dran) of interest, then randomly selecting from the input parameter ranges (tau, or tau_1/tau_2, or tau_kww/beta), convolving with the measured IRF loaded in from the IRF ASCII file, and comparing the RMSE of the convolution of the IRF with the simulated decay using the choice of fit parameters to the experimental lifetime decay trace loaded in by the DataFile. 

Downhill Simplex is used for optimization of the fit function from this point forward to shrink the range of input parameters with each 60-second run such that the maximum values are those that minimized the RMSE over as many Monte Carlo random samples of input parameters can be done in 60 seconds (dependent on local machine). After three 60-second intervals of Monte Carlo least-squares fitting, the fit either converges to the best fit parameters for the model, or exits with the current parameters.

[Least Squares Fitting code](https://github.com/lcgroff2/matlab/tree/master/picofit)
