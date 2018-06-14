
% Batch (background) version of fit routine.

% to run, type: nohup octave fitbatch.m >& run.log &

% User must edit the following lines.

% select data range:
%dran=250:1550;

global IRFx IRFy ;

DataFile='0727-10.asc';
IRFFile='0727-11.asc';

dt=6.7;                 % dwell time in ps (picoseconds per bin)
fixtaus=[10];    % fixed time constants
guessT0=5;    % Guess for time zero
guesstaus=[3100];   %guesses for time constants
global MAXTIME
MAXTIME=10;   %Length of time for Monte Carlo fit.
SaveFile = [DataFile '.fit.mat']; % where to save fit results.
dran=2000:4000;
DoBoxcar=1;
DoErrorbar=0;

% end of user-configurable section

% load datafiles
DataY=loadpico(DataFile);
DataX=(1:length(DataY))*dt;

IRFy=loadpico(IRFFile);
IRFx=DataX;


plot(DataX(dran),DataY(dran)./max(DataY(dran)),IRFx(dran),IRFy(dran)./max(IRFy(dran)))


disp('If plot looks good, press return');
pause

DataX=DataX(dran);
DataY=DataY(dran);
IRFx=IRFx(dran);
IRFy=IRFy(dran);

%set initial T0 to max of IRF
[maxval,maxidx]=max(IRFy);
DataX=DataX-DataX(maxidx);
IRFx=DataX;

if DoBoxcar==1,
  [IRFx,IRFy]=myboxcar(IRFx,IRFy,16);
  [DataX,DataY]=myboxcar(DataX,DataY,16);
end

fprintf(1,'using random number fitting for %4.2f min.\n',MAXTIME/60);

disp('fitting with 25% boundaries on TAU, 20 ps bound on T0')


[pp,yfit,amp]=picofit(DataX,DataY,[(-10+guessT0) (10+guessT0)],fixtaus,0.75*guesstaus,1.25*guesstaus) ;


disp('initial Monte Carlo fit results')
ppmc=pp

chisq=sum((DataY-yfit*amp).^2)/length(DataY)

mcchisq=chisq;

disp('refining fit result using 5 percent boundaries')

guesstaus=pp(2:length(pp));
guessT0=pp(1);

[pp,yfit,amp]=picofit(DataX,DataY,[(-5+guessT0) (5+guessT0)],fixtaus,0.95*guesstaus,1.05*guesstaus) ;
ppsim=pp

chisq=sum((DataY-yfit*amp).^2)/length(DataY)

disp('press any key')

disp('refining fit result using 2 percent boundaries')

guesstaus=pp(2:length(pp));
guessT0=pp(1);

[pp,yfit,amp]=picofit(DataX,DataY,[(-2+guessT0) (2+guessT0)],fixtaus,0.98*guesstaus,1.02*guesstaus) ;
ppsim=pp

disp('refining fit result using 1 percent boundaries')

guesstaus=pp(2:length(pp));
guessT0=pp(1);

[pp,yfit,amp]=picofit(DataX,DataY,[(-1+guessT0) (1+guessT0)],fixtaus,0.99*guesstaus,1.01*guesstaus) ;
ppsim=pp



disp(['T0 is ' num2str(pp(1))])
disp('time constants are')
pp(2:length(pp))

bestfittaus=pp(2:length(pp));
bestfitT0=pp(1);
bestyfit=yfit;
bestamp=amp;



if DoErrorbar==0,
  return
end

disp('now doing error bar estimation')

numfakes=20;

theodata=yfit*amp;

paramstore=zeros([length(pp) numfakes]);

for cnt1=1:numfakes,

  FakeY=shotnoise(theodata);

  % generate fake data with poisson noise
  % note: this is a bad noise generator, **need to fix**
  


  [pp,yfit,amp]=picofit(DataX,FakeY,[(bestfitT0-2) (bestfitT0+2)],fixtaus,0.9*bestfittaus,1.1*bestfittaus);

  pp
  paramstore(:,cnt1)=pp(:);

end


save(SaveFile)

disp('final T0 is ')
bestfitT0
disp('final time constants are ')
bestfittaus

disp('error bars for T0, Tau1, Tau2...  are');
std(paramstore')

plot(DataX,DataY,'o',DataX,theodata,DataX,DataY-theodata)

