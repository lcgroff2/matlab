%timebandwidth.m - Louis Groff 12/12/2016: Calculate time-bandwidth product
%for 50 and 70 fs pulses to determine transform-limited bandwidth of OPA at
%a given center frequency.

wav = input('DFG OPA center wavelength in nm: ');
cm1 = (1/wav)*1e7;
Hztocm1 = 3.33e-11; %1 Hz = 3.33e-11 cm-1

%Gaussian pulse shape:
dEdt_gauss = 0.44; %transform-limited time-bandwidth product of gaussian shaped pulse
dt = 30:2:80;
dt = dt.*1e-15;
dE = dEdt_gauss./dt;
dEcm1 = dE.*Hztocm1;
low = cm1-dEcm1./2;
high = cm1+dEcm1./2;
figure
subplot(2,1,1);
plot(dt,dEcm1,'ok')
title('Gaussian Transform-Limited Bandwidth for given Pulsewidths')
xlabel('pulsewidth (fs)')
ylabel('bandwidth (cm^-^1)')
axis tight
subplot(2,1,2);
plot(dt,low,'or',dt,high,'ob')
title('Gaussian Half-Max Frequencies for given Pulsewidths')
xlabel('pulsewidth (fs)')
ylabel('frequency (cm^-^1)')
axis tight

%Print results:
fprintf(1,'\n')
fprintf(1,'******************************************\n')
fprintf(1, 'For a DFG output wavelength of %4i nm \n',wav)
fprintf(1, 'The center frequency is: %4.0f cm-1\n',cm1)
fprintf(1,'******************************************\n')
