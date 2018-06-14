
% fadfitrun--automated picosecond fit routine.  Changed for anisotropy decay.
% estimates the exponential and kww algebraic amplitudes, instead of "integrated relative amplitudes".
% Also skips "refining" steps, so you should pick fairly narrow search ranges or it will not give good results.

% You should test whether your range is narrow enough and MAXTIME long enough by running the script 
% 2-3 times on the same set of data and comparing the fit results.

% JDM 5/2013

% based on kwwfitrun, JDM 4/2008


global IRFx IRFy MAXTIME;

% User must edit the following lines.

DataFile='jyf071613-n2-para-s1.asc'; % don't get the data and irf files
                            % mixed up, or strange results occur.
IRFFile='jyf071613-n2-mag-i1.asc'; % to get more fit results and possibly
                           % better statistics, try using two
                           % different IRF scans (the one before
                           % the sample scan and the one after)

dt=2.6333;                 % dwell time in ps (picoseconds per
                          % bin). This should be re-measured
                          % periodically.

SaveFile = [DataFile '.fit.mat']; % where to save fit results.
dran=1600:1800; % Data range.  Carefully edit to include 100-300 ps
                % before t0, and to not go much past 4*tau.  Try
                % adjusting dran to test fit: if the fit is robust,
                % changing dran by 50 bins or so on either end
                % should only change the t
                % time constant by a few 
                % percent or less.
ReverseMode=1; % set to "1" for reverse-mode data (PIN sent to START input on TAC)

Rebinlvl=0; % for re-binning data (boxcar, or downsampling), esp. useful for
            % increasing the SNR for samples with long lifetimes.
            % Set Rebinlvl=0 for no re-binning.  Set to 5 for x5
            % reduction in the number of data points.
DoErrorbar=0;  % not implemented yet.  For now, get your error bars
               % by fitting several different sample and IRF runs.
clear fitinfo
curfun=0;

% the "t0 shift" section -- sometimes the IRF needs to be shifted
% by several ps.  Improper t0 shift leads to a sharp feature in the
% residuals that looks like the derivative of the IRF
fitinfo.t0lbound=-70;
fitinfo.t0ubound=50;
fitinfo.t0curval=-30;

% exponential section -- uncomment following block of code to fit an
% exponential, and edit the bounds and curval (initial guess)
curfun=curfun+1;
fitinfo.fun{curfun}.function = 'exp'; % example of exponential
fitinfo.fun{curfun}.lbound=[5];
fitinfo.fun{curfun}.ubound=[2000];
fitinfo.fun{curfun}.curval=[300];

% second exponential (bi-exponential) section -- 
% uncomment following block of code add fit a second
% exponential, and edit the bounds and curval (initial guess)
% BE CAREFUL THAT THE BOUNDS AND GUESS DO NOT OVERLAP WITH THOSE OF
% THE FIRST EXPONENTIAL, OTHERWISE DIVISION BY ZERO CAN RESULT!!!
% curfun=curfun+1;
% fitinfo.fun{curfun}.function = 'exp'; % example of exponential
% fitinfo.fun{curfun}.lbound=[3107.92817];
% fitinfo.fun{curfun}.ubound=[3107.92817];
% fitinfo.fun{curfun}.curval=[3107.92817];

% kww section -- uncomment following block of code to fit a
% kww, and edit the bounds and curval (initial guess)
curfun=curfun+1;
fitinfo.fun{curfun}.function = 'kww'; % can be 'exp', 'kww', 'offset',
                                 % 'linear', 'hside', or 'irf'
fitinfo.fun{curfun}.lbound=[4.33333 0.2998]; % first value is the tau
                                      % parameter, second is beta
fitinfo.fun{curfun}.ubound=[4.33333 0.2998]; % first value is the tau
                                        % parameter, second is beta
fitinfo.fun{curfun}.curval=[4.33333 0.2998]; % first value is the tau
                                      % parameter, second is beta


% offset section -- uncomment following block of code to add a
% constant offset.  This is sometimes necessary due to "dark
% counts", or maybe due to overlap from the previous pulse.
% Has "dummy" bounds and curval
% curfun=curfun+1;
% fitinfo.fun{curfun}.function = 'offset';
% fitinfo.fun{curfun}.lbound=0;
% fitinfo.fun{curfun}.ubound=0;
% fitinfo.fun{curfun}.curval=0;

% heaviside section -- uncomment following block of code to add a
% heaviside, which acts like an exponential with a very long time
% constant.
% Has "dummy" bounds and curval
%curfun=curfun+1;
%fitinfo.fun{curfun}.function = 'hside';
%fitinfo.fun{curfun}.lbound=0;
%fitinfo.fun{curfun}.ubound=0;
%fitinfo.fun{curfun}.curval=0;

% irf section -- uncomment following block of code to add an
% irf, which acts like an exponential with a very short time
% constant.  Also useful to model the effect of scattered laser light somehow getting past the filters.
% Has "dummy" bounds and curval
%curfun=curfun+1;
%fitinfo.fun{curfun}.function = 'irf';
%fitinfo.fun{curfun}.lbound=0;
%fitinfo.fun{curfun}.ubound=0;
%fitinfo.fun{curfun}.curval=0;

MAXTIME=30;   %Length of time for Monte Carlo fit.
% end of user-configurable section

% load datafiles
DataY=loadpico(DataFile);
DataX=(1:length(DataY))*dt;

IRFy=loadpico(IRFFile);
IRFx=DataX;

if ReverseMode==1
    DataY=DataY(end:-1:1); % flip it
    IRFy=IRFy(end:-1:1); % flip it
end

figure(1);
plot(DataX(dran),DataY(dran)./max(DataY(dran)),IRFx(dran), ...
     IRFy(dran)./max(IRFy(dran)))
axis tight

disp('If plot looks good, press return, else hit ctrl-c and change dran');
pause

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

% Some sanity checking on the inputs
if maxval<10000
  disp('***WARNING*** IRF peak value is small (less than 10,000).');
  disp('check "dran"');
end

if DataX(1)<-400
  disp('Warning: a lot of data before t0 is included');
  disp('You might want to trim "dran"');
end

if Rebinlvl>1
  disp('Warning: "Rebinlvl" set.  If lifetime < 1 ns, you might not')
  disp('want to rebin data.  To disable rebinning data, set Rebinlvl')
  disp('to 0 or 1.');
end

  
if Rebinlvl>1,
  [IRFx,IRFy]=myboxcar(IRFx,IRFy,Rebinlvl); % the third parameter is the
                                            % bin size
  [DataX,DataY]=myboxcar(DataX,DataY,Rebinlvl);
end

if length(DataX)>950
  fprintf(1,'*** Warning *** Large Data set (>950 points).\n');
  fprintf(1,'Fitting will converge slowly.  Consider either\n')
  fprintf(1,'increasing MAXTIME, narrowing "dran",\n');
  fprintf(1,'or increasing Rebinlvl\n');
end

  
fprintf(1,'using random number fitting for %i sec.\n',MAXTIME);

[retfitinfo,yfit,amp]=kwwfit(DataX,DataY,fitinfo) ;
fitinfo=retfitinfo;
disp('initial Monte Carlo fit results:')

yfitsum=yfit*amp;
yfitsum=yfitsum(:);

fprintf(1,'RMS residual %.5g percent\n',sqrt(sum((DataY(:)-yfitsum).^2))/length(DataY)/mean(DataY)*100);
fprintf(1,'T0 shift is %5.3f \n',fitinfo.t0curval);
for cnt2=1:length(fitinfo.fun),
  fprintf(1,'function:  %s,  parameters:  ',fitinfo.fun{cnt2}.function);
  fprintf(1, '%8.5f ',fitinfo.fun{cnt2}.curval);
  fprintf(1,'\n');
end

fitinfo0=fitinfo;

disp('refining using jdmmin3 simplex-ish minimization, step 1...')
jdmmin3


% calculate new yfit, amp, etc
[sqerror,yfit,amp]=kwwfit_fun(fitinfo);
yfitsum=yfit*amp;
yfitsum=yfitsum(:);
%figure(1)
%plot(DataX,DataY,'o',DataX,yfit*amp,'g',DataX,yfit*amp-DataY(:))


fprintf(1,'\n======================================\n');
fprintf(1,'  FINAL FIT VALUES\n');
fprintf(1,'======================================\n');
fprintf(1,'\n');
fprintf('dt is %6.4f ps per bin (raw data)\n',dt);
if Rebinlvl>1,
  fprintf(1,'re-binning, with bin width of %i bins\n',Rebinlvl);
  fprintf(1,'dt, taking into account re-binning, is %6.3f ps per bin\n',dt*Rebinlvl);
end
fprintf(1,'number of data points used in fitting: %i\n',length(DataX));
fprintf(1,'RMS residual %.5g percent\n',sqrt(sum((DataY(:)-yfitsum).^2))/length(DataY)/mean(DataY)*100);
fprintf(1,'T0 shift is %5.3f \n',fitinfo.t0curval);
for cnt=1:length(fitinfo.fun),
  fprintf(1,'function:  %s,  parameters:  ',fitinfo.fun{cnt}.function);
  fprintf(1, '%6.5f ',fitinfo.fun{cnt}.curval);
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

disp('Estimated raw amplitudes of each function:')

disp('(a negative amplitude associated with an exponential')
disp('likely indicates an exponential rise)')

for cnt=1:length(fitinfo.fun)
    if fitinfo.fun{cnt}.function == 'kww';
        kwwfun=exp(-(DataX/fitinfo.fun{cnt}.curval(1)).^fitinfo.fun{cnt}.curval(2));
        kwwfun(find(DataX<0))=0; % kww starts at t=0
        kwwint=sum(kwwfun); % times dt for real int, but we skip that here because it cancels
        fitint=sum(bestfitcurves(:,cnt));
        kwwamp=fitint/kwwint; % by algebra and calculus rules, the ratio of the integrals should give the amplitude 
        fprintf(1,'KWW raw amplitude #%i: %.6g counts\n',cnt,kwwamp);
    end
    if fitinfo.fun{cnt}.function == 'exp';
        expfun=exp(-DataX/fitinfo.fun{cnt}.curval(1));
        expfun(find(DataX<0))=0; % exp starts at t=0
        expint=sum(expfun); % times dt for real int, but we skip that here because it cancels
        fitint=sum(bestfitcurves(:,cnt));
        expamp=fitint/expint; % by algebra and calculus rules, the ratio of the integrals should give the amplitude 
        fprintf(1,'EXP raw amplitude #%i: %.6g counts\n',cnt,expamp);
    end
end
        

% clean out the workspace a bit
clear cnt cnt2 cnt3 funran DoBoxcar DoErrorbar
clear maxidx maxval mcchisq retfitinfo t0ran tmplbound
clear tmpubound numrefsteps stepfactor yfitsum
