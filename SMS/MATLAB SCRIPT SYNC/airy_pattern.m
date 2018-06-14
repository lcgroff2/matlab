x = (-1:0.01:1);
lambda = 630; %nm
NA = 1.25;

%function airy=airy_pattern(x,lambda,NA)
lambda=lambda/1000; % lambda in nm, X in microns
k=2*pi/lambda;
idx=find(x==0);
x(idx)=1e-4;
airy=(2 * besselj(1,k * x * NA) ./ (k * x * NA) ).^2;
airy=real(airy);

plot(x,airy)