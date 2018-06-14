FWHM = 290;
a = 65.85;
b = 101;

sCMOS_quanteff = 0.57; % This is the quantum efficiency of the sCMOS camera
collect_eff = 0.03; % This is the collecting efficiency of the microscope 
conv_gain = 0.6; % Be sure to check ADC setup & camera datasheet for conversion gain

sigma = FWHM/(2*sqrt(2*log(2)));
fx = @(x) sqrt(sigma^2/x+a^2/12/x+8*pi*sigma^4*b^2/a^2/x^2);


fplot(fx,[0.1,3])
view(90, -90)
ylabel('Number of photon')
xlabel('Tracking uncertainty (nm)')
grid on

LineH = get(gca, 'Children');
x = get(LineH, 'XData');
y = get(LineH, 'YData');

xaty05ph = interp1(x, y, 0.5);
xaty1ph = interp1(x, y, 1);
xaty2ph = interp1(x, y, 2);

xaty05ct = xaty05ph / (conv_gain/sCMOS_quanteff/collect_eff);
xaty1ct = xaty1ph / (conv_gain/sCMOS_quanteff/collect_eff);
xaty2ct = xaty2ph / (conv_gain/sCMOS_quanteff/collect_eff);

delta = texlabel('delta');

title(sprintf('\nTracking uncertainty, %s\nFluorescence spot size (FWHM) = %3.0f nm.\n%s(%1.2eph=%1.2ecnt)=2  %s(%1.2eph=%1.2ecnt)=1  %s(%1.2eph=%1.2ecnt)=0.5\n', delta,FWHM,delta,xaty2ph,xaty2ct,delta,xaty1ph,xaty1ct,delta,xaty05ph,xaty05ct))

