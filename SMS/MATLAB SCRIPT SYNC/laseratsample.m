
% laseratsample.m -- estimate the intensity profile at the sample plane

% INPUT
LaserPower = 10; % laser power, in microwatts, at back of microscope
PixelSize = 65.85; % pixel size, in nm (might not use this parameter)
%ObjectiveTransmission = 0.1; % estimated/typical 100X 1.25NA objective transmission
SpotFWHM = 36*PixelSize; %5000; % laser spot size, in nm

%%
Sig = SpotFWHM/(2*sqrt(2*log(2)));
SamplePower = LaserPower*1e-6; %ObjectiveTransmission*1e-3; % power at sample plane in W

% for plotting
dx = 10; % dx, in nm
X = -10000:dx:10000; % plotting range
[XX,YY] = meshgrid(X,X);

Z = 1/(2*pi*Sig^2)*exp(-(XX.*XX+YY.*YY)/2/Sig/Sig);
%Z = ng2d(Sig,XX,YY);  % since Sig has units of nm, Z has units of nm^(-2).
%disp('checking to see if ng2d.m is a properly normalized 2D gaussian:')
% should be 1.00 if ng2d is properly normalized:
if (round(sum(sum(Z))*dx*dx*1e6))/1e6 ~= 1
    disp('ng2d.m is not properly normalized in 2D gaussian.')
    break
end

Isample = SamplePower*Z*1e14; % Isample is in W/cm^2
imagesc(X,X,Isample);colorbar
title(sprintf('\nLaser intensity (W/cm^2)\n\nLaser power = %i uW. Laser spot size (FWHM) = %4.0f nm.\nEstimated intensity peak = %3.0f W/cm^2\n\n', LaserPower, SpotFWHM,max(max(Isample))))
%title('Intensity, W/cm^2')
axis square
xlim([-5000 5000])
ylim([-5000 5000])
xlabel('nm')
ylabel('nm')

%fprintf(1,'For a laser power of % 4.2e mW, a laser spot size (FWHM) of % 4.2e nm and an objective transmission of % 4.2e,\n', LaserPower, SpotFWHM, ObjectiveTransmission);
%fprintf(1,'the estimated peak intensity is % 4.2e W per square cm\n',max(max(Isample)));