%Sigmoid simulation and Van't Hoff analysis based on IGOR Pro parameters (LCG):

x = 0:1:100; %Arbitrary temperature vector (�C)

%From IGOR sigmoid fit to real data:
base = 0.11897; %Base of the sigmoid
xhalf = 37.139; %Transition point of the sigmoid
maximum = -0.029159; %Max of the sigmoid
rate = 0.46214; %Decay/Growth rate

%Simulate the sigmoid fit curve based on its parameters:
sig = base+(maximum./(1+exp((xhalf-x)/rate)));
figure
plot(x,sig)
xlabel('Temperature (�C)')
ylabel('Intensity (au)')

%Calculate Keq for Van't Hoff analysis:
n = round(xhalf); %take the midpoint of the curve as an index
if sig(n)-sig(n-1) > 0;
    % For Postive sloping Sigmoids calculate Keq using absolute max and min:
    Keq = (sig-min(sig))./(max(sig)-min(sig));
    figure
    plot(x,Keq)
    xlabel('Temperature (�C)')
    ylabel('Keq')
    axis([x(1) x(end) -0.05 1.05])
else
    %For Negative Sloping Sigmoids calculate Keq taking the absolute max as the
    %"pre-transition" minimum:
    Keq = (sig-max(sig))./(min(sig)-max(sig));
    figure
    plot(x,Keq)
    xlabel('Temperature (�C)')
    ylabel('Keq')
    axis([x(1) x(end) -0.05 1.05])
end

%Generate the Van't Hoff plot:
T = x+273.15; %Need temperature in Kelvin
figure
plot(T.^(-1),log(Keq))
xlabel('1/T (K^-^1)')
ylabel('ln(Keq)')
title('Vant Hoff Plot')