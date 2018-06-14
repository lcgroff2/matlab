function [xout,yout]=boxcar(x,y,n)
%BOXCAR   Reduce the number of points in an X,Y data pair by averaging.
%   [XB,YB]=BOXCAR(X,Y,N) will return XB,YB containing a factor of N fewer
%   points than [X,Y].
x=x(:);
y=y(:);
sz=length(y);
yout=zeros(size(y(1:n:(sz-n+1))));
xout=yout;
for cnt=1:n,
	yout=yout+y(cnt:n:(sz+cnt-n));
   xout=xout+x(cnt:n:(sz+cnt-n));
end
yout=yout/n;
xout=xout/n;

