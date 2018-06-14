% tsta2c2d.m

tic;ac2d=acorr2d(simcmos,100);toc

ac2d=ac2d-ac2d(100,100);
ac2d=ac2d-ac2d(1,1);
imagesc(ac2d);
