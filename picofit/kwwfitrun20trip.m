
% kwwfitrun--automated picosecond fit routine.  JDM 4/2008

% Note: see text at bottom of script for calculating weighted
% time constants.

% User must edit the following lines.

global IRFx IRFy ;

IRFFile='crs102909-3.asc';
DataFile='crs102909-4.asc';


dt=2.65;                 % dwell time in ps (picoseconds per bin)
clear fitinfo
fitinfo.t0lbound=-20;	%-49
fitinfo.t0ubound=80;	%-20
fitinfo.t0curval=37.5;
fitinfo.fun{1}.function = 'kww'; % can be 'exp', 'kww', 'offset', 'linear', 'hside', or 'irf'
fitinfo.fun{1}.lbound=[15 .2];	%[12 .2]
fitinfo.fun{1}.ubound=[100 1];	%[100 1]
fitinfo.fun{1}.curval=[32 .3];
% fitinfo.fun{1}.function = 'exp'; % example of exponential
% fitinfo.fun{1}.lbound=[200];
% fitinfo.fun{1}.ubound=[500];
% fitinfo.fun{1}.curval=[363];
% fitinfo.fun{2}.function = 'exp'; % Add this one if you need to use a second exponential
% fitinfo.fun{2}.lbound=[1];
% fitinfo.fun{2}.ubound=[200];
% fitinfo.fun{2}.curval=[89];
% fitinfo.fun{3}.function = 'exp'; % Add this one if you need to use a third exponential
% fitinfo.fun{3}.lbound=[5];
% fitinfo.fun{3}.ubound=[200];
% fitinfo.fun{3}.curval=[16];

%fitinfo.fun{2}.function = 'offset'; % offset has no nonlinear
                                    % parameters, but need a dummy
                                    % lbound,ubound,curval
%fitinfo.fun{2}.lbound=0;
%fitinfo.fun{2}.ubound=0;
%fitinfo.fun{2}.curval=0;
%fitinfo.fun{3}.function = 'hside'; % heaviside has no parameters
%fitinfo.fun{3}.lbound=0;
%fitinfo.fun{3}.ubound=0;f
%fitinfo.fun{3}.curval=0;
%fitinfo.fun{4}.function = 'irf'; % irf has no parameters
%fitinfo.fun{4}.lbound=0;
%fitinfo.fun{4}.ubound=0;
%fitinfo.fun{4}.curval=0;

global MAXTIME
MAXTIME=45;   %Length of time for Monte Carlo fit.
SaveFile = [DataFile '.fit.mat']; % where to save fit results.
dran=4450:5850;
DoBoxcar=1;
DoErrorbar=0;

% end of user-configurable section

% load datafiles
DataY=loadpico(DataFile);
DataX=(1:length(DataY))*dt;

IRFy=loadpico(IRFFile);
IRFx=DataX;


plot(DataX(dran),DataY(dran)./max(DataY(dran)),IRFx(dran),IRFy(dran)./max(IRFy(dran)))


disp('If plot looks good, press return, else hit ctrl-c and change dran');
% pause

DataX=DataX(dran);
DataY=DataY(dran);
IRFx=IRFx(dran);
IRFy=IRFy(dran);

% subtract off the background as the average of first 2 data points
DataY=DataY-mean(DataY(1:2));
IRFy=IRFy-mean(IRFy(1:2));

%offset X axis to max of IRF
[maxval,maxidx]=max(IRFy);
DataX=DataX-DataX(maxidx);
IRFx=DataX;

if DoBoxcar==1,
  [IRFx,IRFy]=myboxcar(IRFx,IRFy,2);
  [DataX,DataY]=myboxcar(DataX,DataY,2);
end

fprintf(1,'using random number fitting for %4.2f min.\n',MAXTIME/60);


[retfitinfo,yfit,amp]=kwwfit(DataX,DataY,fitinfo) ;



disp('initial Monte Carlo fit results:')
retfitinfo
retfitinfo.fun{:}

yfitsum=yfit*amp;
yfitsum=yfitsum(:);
chisq=sum((DataY(:)-yfitsum).^2)/length(DataY)

mcchisq=chisq;

disp('refining...')

stepfactor = 1/4; % decrease range by this factor for each
                  % refinement step
numrefsteps = 2;     % number of refinement steps

for cnt=1:numrefsteps
  fitinfo=retfitinfo; % use last best guess as starting point and
                      % calculate new narrower bounds
  t0ran=fitinfo.t0ubound-fitinfo.t0lbound;
  t0ran=t0ran*stepfactor;
  fitinfo.t0lbound=fitinfo.t0curval-t0ran/2;
  fitinfo.t0ubound=fitinfo.t0curval+t0ran/2;
  for cnt2=1:length(fitinfo.fun)
    funran=fitinfo.fun{cnt2}.ubound-fitinfo.fun{cnt2}.lbound;
    funran=funran*stepfactor;
    for cnt3=1:length(funran),
      tmplbound=fitinfo.fun{cnt2}.curval(cnt3)-funran(cnt3)/2;
      if tmplbound>fitinfo.fun{cnt2}.lbound(cnt3),
		fitinfo.fun{cnt2}.lbound(cnt3)=tmplbound;  % don't go lower
                                                   % than initial
                                                   % lower bound
      end
      tmpubound=fitinfo.fun{cnt2}.curval(cnt3)+funran(cnt3)/2;
      fitinfo.fun{cnt2}.ubound(cnt3)=tmpubound;
    end
  end
  [retfitinfo,yfit,amp]=kwwfit(DataX,DataY,fitinfo) ; % calculate a
                                                      % new fit
  yfitsum=yfit*amp;
  yfitsum=yfitsum(:);
  chisq=sum((DataY(:)-yfitsum).^2)/length(DataY)
  retfitinfo
  retfitinfo.fun{:}
end % end of refinement loop

fitinfo=retfitinfo;



fprintf(1,'T0 is %5.1f \n',fitinfo.t0curval);
for cnt=1:length(fitinfo.fun),
  fprintf(1,'function:  %s,  parameters:  ',fitinfo.fun{cnt}.function);
  fprintf(1, '%6.3f',fitinfo.fun{cnt}.curval);
  fprintf(1,'\n');
end




bestyfit=yfit;
bestamp=amp;
bestfitcurves=zeros(size(yfit));
theodata=yfit*amp;
for cnt=1:length(amp),
  bestfitcurves(:,cnt)=yfit(:,cnt)*amp(cnt);
end



if DoErrorbar==1,

  disp('have not implemented error bar estimation yet')

end

save(SaveFile)

figure(1)
plot(DataX,DataY,'o',DataX,theodata,DataX,DataY(:)-theodata(:))
figure(2)
plot(DataX,bestfitcurves)
title('Individual fit curves')
figure(1)

disp('Integrated relative amplitudes of each function:')

disp('(a negative amplitude associated with an exponential')
disp('likely indicates an exponential rise)')

intamps=sum(bestfitcurves);

relintamps=intamps/sum(intamps)

% for the case of 4 time constants, assuming the first one not of interest
% you would calculate the photon-weighted average time constant as:

% sum(intamps(2:4)/sum(intamps(2:4)).*bestfittaus(2:4))

% clean out the workspace a bit
clear cnt cnt2 cnt3 funran dt DoBoxcar DoErrorbar
clear maxidx maxval mcchisq retfitinfo t0ran tmplbound
clear tmpubound numrefsteps stepfactor yfitsum
