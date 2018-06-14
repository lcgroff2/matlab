function yout=shotnoise(x)

% SHOTNOISE -- an attempt at an efficient shot noise (counting noise)
%              generator.  The True shot noise is most closely
%              approximated using poisson_rnd(), but poisson_rnd()
%              is slow at larger values of "x".  Here we use a
%              gaussian random number generator to approximate
%              the poisson distribution.  For large <n>, the poisson
%              distribution is approximated by a gaussian centered at
%              <n> with a width of <n>^(1/2).  For <n> larger than 20,
%              the poisson distribution and the gaussian distribution
%              are very similar.  Since most of our experiments involve
%              counts above 20, this is a good approximation for this
%              case.

% For example, let's look at the speed of generating 1000 points
% of shot noise, with <n>=400:

% x=ones([1 1000]) .* 400;

% >> tic; y=shotnoise(x); toc
% Elapsed time is 0.000276 seconds.

% >> tic; y=poisson_rnd(x); toc
% Elapsed time is 0.325704 seconds.

% this shows shotnoise() is much faster than poisson_rnd().
yout = x + sqrt(x) .* randn(size(x));

