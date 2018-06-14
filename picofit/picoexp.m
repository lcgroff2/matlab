function yval=simexp(xdat,pp)
  %function yval=simexp(xdat,[b1,b2, ...])
  %
  %simexp: "multi exponential "
  %
  %	Returns a matrix containing multiple exponentials.
  %
  %     example 1: x=0:100;plot(x,simexp(x,[10 20 30 40]))
  %
  %     This example plots 3 exponentials starting at x=10, with
  %     decay constants of 20, 30, and 40.
  %
  %     example 2: plot(x,[1 2 3]*simexp(x,[10 20 30 40])
  %
  %     As above, except plots the sum of the first exponential,
  %     two times the second exponential, and three times the third
  %     exponential.

  %calculate exponential
  xdat=xdat(:)';
  xdat=xdat-min(xdat);
  nparms=length(pp);
  yval=zeros([(nparms) length(xdat)]);
  for cnt=1:nparms,
    yval(cnt,:) = exp(-1/pp(cnt)*xdat);
    %force exponential to zero for negative x
 end
 
