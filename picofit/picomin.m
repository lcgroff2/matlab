
function pout=mcminbnd(F,lbound,ubound)

global MAXTIME
% user set parameters

if exist('MAXTIME')==0,
  MAXTIME=30;
end

pin=(lbound+ubound)/2;

%if length(abs(ubound-lbound)<minbndsz)>0,
%  bdidx=find(abs(ubound-lbound)<minbndsz);
%  ubound(bdidx)=ubound(bdidx)+minbndsz;
%  lbound(bdidx)=lbound(bdidx)-minbndsz;
%end

tic;

bestguess=pin;
bestminval=feval(F,pin);

brange=ubound-lbound;

numtries=0;
while toc<MAXTIME,
  numtries=numtries+1;
  pguess=lbound+(ubound-lbound).*rand(size(lbound));
  curminval=feval(F,pguess);
  if curminval<bestminval,
    bestguess=pguess;
    bestminval=curminval;
    % the following should be commented out for strict bounded search
    %lbound=pguess-brange/2;
    %ubound=pguess+brange/2;
  end
end
numtries
pout=bestguess;

