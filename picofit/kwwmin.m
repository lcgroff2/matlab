
function bestguess=kwwmin(fitinfo)

global MAXTIME
% user set parameters


if exist('MAXTIME')==0,
  MAXTIME=30;
end

%pin=(lbound+ubound)/2;

%if length(abs(ubound-lbound)<minbndsz)>0,
%  bdidx=find(abs(ubound-lbound)<minbndsz);
%  ubound(bdidx)=ubound(bdidx)+minbndsz;
%  lbound(bdidx)=lbound(bdidx)-minbndsz;
%end

tic;

bestguess=fitinfo;
fitguess=fitinfo;
bestsqerr=kwwfit_fun(fitinfo);

%brange=ubound-lbound;

numtries=0;
while toc<MAXTIME,
  numtries=numtries+1;
  % first generate a guess for function parameters
  for cnt=1:length(fitguess.fun)
    fitguess.fun{cnt}.curval = fitguess.fun{cnt}.lbound+(fitguess.fun{cnt}.ubound-fitguess.fun{cnt}.lbound).*rand(size(fitguess.fun{cnt}.lbound));
  end
  % now generate a guess for t0
  fitguess.t0curval=fitguess.t0lbound+(fitguess.t0ubound-fitguess.t0lbound)*rand(1);
  cursqerr=kwwfit_fun(fitguess);
  if cursqerr<bestsqerr,
    bestguess=fitguess;
    bestsqerr=cursqerr;
    % the following should be commented out for strict bounded search
    %lbound=pguess-brange/2;
    %ubound=pguess+brange/2;
  end
end
fprintf('number of function evaluations %i\n',numtries);


