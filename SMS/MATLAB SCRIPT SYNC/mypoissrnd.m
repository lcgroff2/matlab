
function y=mypoissrnd(lambda)

% simulates poisson random number generator, but is faster for
% larger lambda.  Not accurate for lambda less than about 10.

y=round(lambda+randn(size(lambda)).*sqrt(lambda));
y(find(y<0))=0;
lowidx=find(lambda<10);
y(lowidx)=poissrnd(lambda(lowidx));
