
% pltnr.m -- script to calculate and plot exciton distribution, n(r), from exlatsim.

% can plot either the excitons that undergo energy transfer, or those that don't,
% depending on whether you use etidxstore or decidxstore

nruns=1;
allpos1=[];
allpos2=[];

for runidx=1:nruns;
    exlatsim
    allpos1=cat(2,allpos1,etposstore); % for positions of excitons that undergo ET
    allpos2=cat(2,allpos2,decposstore); % for positions of excitons that undergo decay

end

save rdisttmp.mat

% next, convert each XYZ position into an R value
allpos=decposstore;
allR=sqrt( (allpos(1,:)-qposx(1)).^2 + (allpos(2,:)-qposy(1)).^2 + (allpos(3,:)-qposz(1)).^2 );

% then histogram it. (might need a special homemade histogramming function)
xhist=0:.5:7;
Nhist=histc(allR,xhist);

% then calculate normalization factors, i.e., volume factor (see notes).
cubeR=sqrt( (XX-qposx(1)).^2 + (YY-qposy(1)).^2 + (ZZ-qposz(1)).^2 );

cubeR=cubeR(:); % reshape into vector--will histogram them, so order doesn't matter
Jhist=hist(cubeR,xhist);

save rdistresults.mat