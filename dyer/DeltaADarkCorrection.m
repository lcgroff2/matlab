% Script to calculate artificial single point differential absorbance 
%with and without dark current:

%fake A/D counts

sample = 1000; %raw laser intensity no absorbance
ref = 100;     %reference laser intensity
signal = 500;  %assume 500 count shift in detected intensity due to 
               %absorbance change
dark = 50;    %dark current

%with no dark current correction (assume a bleach):
abs_on1 = -log10((sample+signal+dc)/(ref+dc));
abs_off1 = -log10((sample+dc)/(ref+dc));
delta_A1 = abs_on1 - abs_off1;

%with dark current correction:
abs_on2 = -log10((sample+signal)/ref);
abs_off2 = -log10(sample/ref);
delta_A2 = abs_on2 - abs_off2;

delta_A1
delta_A2
dDeltaA = delta_A1-delta_A2