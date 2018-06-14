% cntrfastttest -- test fast centroid calculator "cntrfast"

rand('twister',sum(100*clock)); % randomize seed for random num generator

ix=1:80;
iy=1:60;

bp1=1;
bp2=20;

%xcnt=40.63;
%ycnt=30.17;
width=1.8; % sigma width.  width=1.8 gives roughly 4.25 pixel fwhm
[XX,YY]=meshgrid(ix,iy);


numits=1000;

% make some variables to store results of simulations in
x3dg=zeros([1 numits]);
y3dg=x3dg;
x2dg=x3dg;
y2dg=x3dg;
xcntrd=x3dg;
ycntrd=x3dg;
xfast=x3dg;
yfast=y3dg;
tcntrd=x3dg; % for storing timing results
t2dg=x3dg;
t3dg=x3dg;
tfast=x3dg;% "

H=mkffilt(XX,bp1,bp2); % only need to make the filter fun once

tic; % start timer

for cnt=1:numits,
    xtrue=39+rand(1)*4;
    ytrue=28+rand(1)*4;
    img=1000*exp(-((XX-xtrue).^2+(YY-ytrue).^2)/2/width^2);
    % The above gives roughly 1e5 total photons for width=1.8
    % add some gaussian poisson-like noise:
    img=img+sqrt(img).*randn(size(img));
    % Now add some readout noise (6 electrons) and an offset
    img=img+400+6*randn(size(img));
    peaks=round([xtrue ytrue]);
    % Now fit using 3D gaussian
    %centerg for old fit style, cntr3dg for new fit routine
    t1=toc;
    centers=cntr3dg(img,peaks,5);
    t2=toc;    
    x3dg(cnt)=centers(1)-xtrue;
    y3dg(cnt)=centers(2)-ytrue;
    t3dg(cnt)=t2-t1;
    numphotons=sum(sum(img));
    % Now fit using 2D gaussian
    % there are 3 variants on the 2D gaussian fit: centerg, centerf, centerf3, in order of "sophistication"
    % centerf3 seems to be the fastest, though it could likely be optimized further
    t1=toc;
    centers=centerf3(img,peaks,5);
    t2=toc;    
    x2dg(cnt)=centers(1)-xtrue;
    y2dg(cnt)=centers(2)-ytrue;
    t2dg(cnt)=t2-t1;
    % now do Grier/Crocker/Weeks cntrd for comparison
    %imgbp=bpass(img,1,11);
    imgf=fpass(img,H);
    t1=toc;
    centers=cntrd(imgf,peaks,11); % note their "sz" is 2*sz+1 as compared to cntr3dg
    t2=toc;
    xcntrd(cnt)=centers(1)-xtrue;
    ycntrd(cnt)=centers(2)-ytrue;
    tcntrd(cnt)=t2-t1;
    % now do "fast" centroid code
    t1=toc;
    centers=cntrfast(imgf,peaks,5); % note their "sz" is 2*sz+1 as compared to cntr3dg
    t2=toc;
    xfast(cnt)=centers(1)-xtrue;
    yfast(cnt)=centers(2)-ytrue;
    tfast(cnt)=t2-t1;
end
fprintf(1,'Fitting statistics, %i runs\n',numits);
fprintf(1,'uncertainty from fitting 3d gaussian (cntr3dg) is: %6.4f pixel\n',(std(x3dg)+std(y3dg))/2);
fprintf(1,'uncertainty from fitting 2d gaussian (centerg/f/f3) is: %6.4f pixel\n',(std(x2dg)+std(y2dg))/2);
fprintf(1,'uncertainty from Grier/Crocker cntrd method is: %6.4f pixel\n',(std(xcntrd)+std(ycntrd))/2);
fprintf(1,'uncertainty from cntrfast centroid method is: %6.4f pixel\n',(std(xfast)+std(yfast))/2);
fprintf(1,'cpu time per frame, cntr3dg: %.4g ms\n',mean(t3dg)*1000);
fprintf(1,'cpu time per frame, centerf: %.4g ms\n',mean(t2dg)*1000);
fprintf(1,'cpu time per frame, cntrd: %.4g ms\n',mean(tcntrd)*1000);
fprintf(1,'cpu time per frame, cntrfast: %.4g ms\n',mean(tfast)*1000);
fprintf(1,'number of photons %.2e\n',round(numphotons));