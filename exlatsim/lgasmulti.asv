
% simple script to run lattice gas model multiple times and record quenching efficiency and dynamics.

% the parameters for the run are specified in exlatsim

% unlike etmulti, it doesn't let you switch parameters between runs, and doesn't do fitting to KWW, etc.

nrun=50;
qearray=zeros([1 nrun]);
betavec=zeros([1 nrun]);

for rcnt=1:nrun
    exlatsim2;
    qearray(rcnt)=qe;
    betavec(rcnt)=pp(1);
    if rcnt==1
        picoav=expop; % save the picosecond dynamics
    else
        picoav=picoav+expop; % average the picosecond dynamics
    end
end

tauave=sum(t.*picoav)/sum(picoav);
fprintf('***\noverall average quenching efficiency for %i runs: %.4f\n',nrun,mean(qearray));
fprintf('picosecond decay in t, picoav\n');
fprintf('average tau = %4.3g ps\n',tauave);
fprintf('average beta = %.4f\n',mean(betavec));
plot(t,picoav(min(find(picoav<(0.05*max(picoav))))),'o')

fitexlat