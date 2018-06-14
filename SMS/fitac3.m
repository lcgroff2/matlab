%ac3=jmxcorrslow(simcmos.*simcmos,simcmos); % corr(I^2,I) 3rd-order autocorrelation

dI=simcmos-mean(simcmos); % maybe this will help

ac3f=jmftxcorr(dI.*dI,dI); % trying fft version of my xcorr code for 3rd order

acsq=nac(dI.*dI);

expfun=@(parm,xdata)(parm(1)*exp(-xdata./parm(2))+parm(3));

% acfit=ac3f; % pick either this line or next
acfit=acsq;

iniguess(1)=acfit(1)-acfit(300); % initial guess for exponential amplitude

iniguess(2)=1/ktot; % initial guess for lifetime parameter

iniguess(3)=acfit(300); % initial guess for offset

ran=1:500;
warning off
[outparm,resid,J,Sigma]=nlinfit(t(ran),acfit(ran),expfun,iniguess);
warning on

fity=expfun(outparm,t(ran));

plot(t(ran),acfit(ran),'o',t(ran),fity)

fprintf('ac3ft: A0: %.5g, tau: %.5g, A1: %.5g\n',outparm(1),outparm(2),outparm(3));
