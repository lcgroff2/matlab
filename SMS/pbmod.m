% photobleaching curve with laser fluctuation correction by xiaoli 08/30/2011%

merge = input('Merge to the previous plots? (y/N) ','s');
if merge == 'y' | merge == 'Y';
    disp('-Yes.')
    hold all
elseif isempty(merge) | merge == 'n' | merge == 'N';
    disp('-No.')
    clc
    clear
    close all
else
    disp('-Invalid answer.')
    return
end

cdpwd = pwd;
set(0, 'DefaultAxesFontSize', 14, 'DefaultAxesFontName', 'Times New Roman')
set(0, 'DefaultTextFontSize', 14, 'DefaultTextFontName', 'Times New Roman')
%set(0,'defaulttextfontsize',20);

cc = [0,0,0.75;
    1,0,0;
    0,0.5,0;
    0.75,0.25,0.75;
    0.55,0.27,0.07];

if ispc, %ISPC--True for the PC (Windows) version of MATLAB.ISPC, returns 1 for PC (Windows) versions of MATLAB and 0 otherwise.
    disp('Choose an SPE file.')
    [fnam, pathnam, filterindex] = uigetfile('*.spe', 'Pick an SPE-file');
    %cd(pathnam);
else
    ls *.spe
    fnam = input('enter SPE file name ','s');
end

disp(fnam)
II = loadspe([pathnam fnam]);
[SPE.dim_i,SPE.dim_j,framenum] = size(II);

%%
%plot the first frame of raw spectrum%
figure(1)
peakII = max(II);
[maxII,maxIIloc] = max(peakII);
subplot(2,2,1), plot(II(:,maxIIloc)); hold all
set(gca,'xlim',[0 1024])
title('The frame at highest intensity')
xlabel('X (Pixel)')
ylabel('Intensity (CCD counts)')

%%
%Remove cosmic rays%
for cnt = 1:framenum
    IId(:,cnt) = despike(II(:,cnt));
end

%%
%ADC OFFSET
%determine the ADC 'offset' per pixel by averaging small number of pixels somewhere to the left of the laser, averaging over all frames.%

%Average from x=x1 to x=x2%
laserpx = input('ADC offset will be determined from Pixel 50 to 100. (Y/n) ','s');
if laserpx == 'n' | laserpx == 'N';
    disp('-No.')
    strtpx = input('Start pixel = ');
    endpx = input('End pixel = ');
elseif isempty(laserpx) | laserpx == 'y' | laserpx == 'Y';
    disp('-Yes.')
    strtpx = 50;
    endpx = 100;
else
    disp('-Invalid answer.')
    cd(cdpwd)
    return
end

ADCoffset = round(mean(mean(IId(strtpx:endpx,:)))); % averaging from startpx to endpx for all the frames%

%subtract the offset from all of the frames%
IIbase=IId-ADCoffset;

subplot(2,2,3), plot(IIbase(:,maxIIloc)); hold all
set(gca,'xlim',[0 1024])
title('The frame at highest intensity after subtracting ADC offset')
xlabel('X (Pixel)')
ylabel('Intensity (CCD counts)')

%%
%LASER INTENSITY
%determine the laser intensity, by summing over a few pixels in each frame (we want to integrate 
%over enough photons to get at least 200:1 signal/noise ratio, so at least 4000-5000 total counts above the offset). 
%Be sure to properly correct for ADC offset.%
%Sum from X3 to x4 for each frame%

laserpx = input('Laser intensity will be determined from Pixel 430 to 465. (Y/n) ','s');
if laserpx == 'n' | laserpx == 'N';
    disp('-No.')
    strtlspx = input('-Start laser pixel = ');
    endlspx = input('-End laser pixel = ');
elseif isempty(laserpx) | laserpx == 'y' | laserpx == 'Y';
    disp('-Yes.')
    strtlspx = 430;
    endlspx = 465;
else
    disp('-Invalid answer.')
    cd(cdpwd)
    return
end

for cnt = 1:framenum
    laserbase(cnt) = sum(IIbase(strtlspx:endlspx,cnt));
end

% divide the array of laser intensity at each frame by its mean, to normalize it to 1, like: norm_laser = cnts_laser/mean(cnts_laser) %
norm_laserbase = laserbase/max(laserbase); 
subplot(2,2,2), plot(norm_laserbase); hold all
title('Normalized laser fluctuation')
xlabel('Framenumber')
ylabel('Normalized laser I')

%%
%FLUORESCENCE INTENSITY
%similarly, pick a few pixels for the fluorescence, and correct for offset, to generate fluorescence photobleach kinetics trace.%

%Sum
flpx = input('Fluorescence intensity will be determined from Pixel 500 to 1024. (Y/n) ','s');
if flpx == 'n' | flpx == 'N';
    disp('-No.')
    strtflpx = input('-Start fluorescence pixel = ');
    endflpx = input('-End fluorescence pixel = ');
elseif isempty(flpx) | flpx == 'y' | flpx == 'Y';
    disp('-Yes.')
    strtflpx = 500;
    endflpx = 1024;
else
    disp('-Invalid answer.')
    cd(cdpwd)
    return
end

for cnt=1:framenum
    fluobase(cnt)=sum(IIbase(strtflpx:endflpx,cnt));
end

subplot(2,2,4), plot(fluobase); hold all
title('Summed fluorescence intensity before laser correction')
xlabel('Framenumber')
ylabel('Summed I')

%divide the photobleach kinetics by the laser intensity and normalize photobleaching curve%

fluocorr = fluobase./norm_laserbase;
fluomaxcheck = max(fluocorr); %If this value is 'inf', there is a value of norm_laserbase = 0

%%
%TIME SCALE
tintvl = input('Integration time interval (sec/frame) = ');
t = 0:tintvl:(framenum-1)*tintvl;
I = fluocorr;

modulation = input('Laser modulation analysis? (Y/n) ','s');
if isempty(modulation) | modulation == 'y' | modulation == 'Y';
    disp('-Yes.')
    figure(2)
    subplot(3,1,1), plot(t,I); hold all
    title('Normalized photobleaching curve')
    xlabel('Time (s)')
    ylabel('Summed I / Laser I (I")')
elseif modulation == 'n' | modulation == 'N';
    disp('-No.')    
    bleaching = input('Merge photobleaching to modulation results? (y/N) ','s');
    if bleaching == 'y' | bleaching == 'Y';
        disp('-Yes.') 
        %Differentiate peak of laser intensity
        diffls = abs([diff(laserbase),0]);
        [peakls,locsls] = findpeaks(diffls,'minpeakheight',1.5*std(diffls)+mean(diffls));
        %Crop pixel
        tminpx = t(I==I(1,locsls(1,1)+2)) / tintvl; %Pixel of tmin at Imax
        tbleach = t(1,tminpx:end) - t(I==I(1,locsls(1,1)+2)) + tintvl; tbleach(1,1) = 0; %Crop pixels from tminpx to tmaxpx. Relative the min to zero.
        Ibleach = I(1,tminpx:end); %Crop pixels of I as the same length of t.
        %Normalize time interval
        inputnmltime = input('Time normalization? ');
        nmltbleach = tbleach*inputnmltime;
        %Plot
        figure(2)
        subplot(3,1,3);
        plot(nmltbleach,Ibleach/I(1,locsls(1,1)+2)); hold all
        cd(cdpwd)
        return    
    elseif isempty(bleaching) | bleaching == 'n' | bleaching == 'N';
        disp('-No.') 
        figure(5)
        plot(t,I/max(I)); hold all
        title('Normalized photobleaching curve')
        xlabel('Time (s)')
        ylabel('Summed I / Laser I')
        cd(cdpwd)
        return    
    else
        disp('-Invalid answer.')
        cd(cdpwd)
        return
    end
else
    disp('-Invalid answer.')
    cd(cdpwd)
    return
end

%%
%MODULATION
period = input('Period (sec) = ');

%Find the first period
nmllaserbase = laserbase/max(laserbase);
diffls = abs([diff(nmllaserbase),0]);
[peakls,locsls] = findpeaks(diffls,'minpeakdistance',period/(2*tintvl)*0.70);
if locsls(1,1) < 50 ;
    locsls = locsls(1,2:end);
end

%Crop pixel (six periods)
tminpx = t(I==I(1,locsls(1,1))) / tintvl; %Pixel of tmin at Imax
tmaxpx = tminpx + 6*period / tintvl - 1; %Pixel of tmax with six periods.
t6p = t(1,round(tminpx:tmaxpx)) - t(I==I(1,locsls(1,1))) + tintvl; t6p(1,1) = 0; %Crop pixels from tminpx to tmaxpx. Relative the min to zero.
I6p = I(1,round(tminpx:tmaxpx)); %Crop pixels of I as the same length of t.

%Max I
%startI = mean(I(1,round(locsls(1,1)-35:locsls(1,1)-5))); %Average from threshold-5 to threshold-35 (N=30)
startI = I(1,locsls(1,1)+2);
thresholdI = I(1,locsls(1,1));
%maxI = max(I(1,round(locsls(1,1)-2/tintvl):1:round(locsls(1,1)+2/tintvl)));
%maxI = t(I==thresholdI)/tintvl+1;

%Laser
ls6p = nmllaserbase(1,round(tminpx:tmaxpx));
%Diff laser
diffls6p = abs([diff(ls6p),0]);

%Plot I and laser v.s. time
subplot(3,1,2);
[ax,Ilsvst] = plotyy(t6p,I6p/startI,t6p,ls6p); hold all %Normalize I.
set(ax(1),'ylim',[0 1.4],'ytick',[0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8],'box','off','ygrid','on'); 
set(ax(2),'ylim',[0 5],'ytick',[0 1]); 
set(ax,'xlim',[0 6*period]);
set(get(ax(1),'ylabel'),'string','Normalized I"');
set(get(ax(2),'ylabel'),'string','Normalized Laser');
xlabel('Time (s)')
title('Normalized photobleaching curve for first six periods, Normalized laser, and Differentiated laser')

%Plot differentiated laser
line(t6p,diffls6p,'color','r'); hold all

%Peak finding
[peak,locs] = findpeaks(diffls6p,'minpeakdistance',period/(2*tintvl)*0.70,'minpeakheight',0.3); %the distance between peak is 70% of the period (before normalized).
for i=2:length(locs)-1;
    tintv(i-1) = locs(1,i+1)-locs(1,i);
end

%Crop five period
startlocs = 1;
t5pminpx = locs(1,startlocs);
t5pmaxpx = locs(1,startlocs)+5*2*mean(tintv)+1;
t5p = t6p(1,ceil(t5pminpx:t5pmaxpx)) - t6p(1,ceil(t5pminpx)); %Crop pixels from tminpx to tmaxpx. Relative the min to zero.
I5p = I6p(1,ceil(t5pminpx:t5pmaxpx)); %Crop pixels of I as the same length of t.

%Normalize time interval
nmlfactor = period/(tintvl*2*mean(tintv));
nmlt5p = t5p*nmlfactor;

%Plot
subplot(3,1,3);
plot(nmlt5p,I5p/startI); hold all
set(gca,'xlim',[0 5*period],'xtick',[0:period/2:5*period],'xgrid','on')
title('Normalized photobleaching curve of period 2 to 6')
xlabel('Normalized Time (s)')
ylabel('Normalized I"')

%Cut the intensity each period
thalfpintv = mean(tintv);
figure(3)
clrr = 1;
for i = 1:1:10;
    Ihalfp(i,:) = I5p(1,ceil((i-1)*thalfpintv+2):ceil((i-1)*thalfpintv+2)+floor(thalfpintv-2));
    cutthalfp(i,:) = nmlt5p(1,ceil((i-1)*thalfpintv+2):ceil((i-1)*thalfpintv+2)+floor(thalfpintv-2));
    if mod(i,2) %Photobleaching 1st half-period
        nmlthalfp(i,:) = cutthalfp(i,:) + (1-i)*period/2;
    else %Reverse 2nd half-period
        nmlthalfp(i,:) = cutthalfp(i,:) + (2-i)*period/2;
    end
    
    subplot(1,2,1)
    plot(nmlthalfp(i,:),Ihalfp(i,:)/startI,'color',cc(clrr,:)); hold on
    title('Normalized photobleaching curve of period 2 to 6')
    xlabel('Normalized Time (s)')
    ylabel('Normalized I"')
    set(gca,'xlim',[0 period],'xtick',[0:period/2:period],'xgrid','on')
    
    for j = 6*tintvl:tintvl/2:period/2-6*tintvl;
        if mod(i,2)
            I5psplit(i,round(j*2/tintvl-11)) = interp1(nmlthalfp(i,:),Ihalfp(i,:),j,'linear');
        else
            I5psplit(i,round(j*2/tintvl-11)) = interp1(nmlthalfp(i,:),Ihalfp(i,:),j+period/2,'linear');
        end
    end
    
    if mod(i,2)
    else 
        clrr = clrr +1;
    end

end

t5psplit1sthalf = 6*tintvl:tintvl/2:period/2-6*tintvl;
t5psplit2ndhalf = t5psplit1sthalf + period/2;
%Sum the intensity of cut half-periods
meanI1sthalf(1,:) = mean(I5psplit(1:2:9,:),1);
meanI2ndhalf(1,:) = mean(I5psplit(2:2:10,:),1);

subplot(1,2,2)
[fitgrph1st,fit1details] = fit(t5psplit1sthalf(1,:)',meanI1sthalf(1,:)'/startI,'exp2')  
fitplot1st = plot(fitgrph1st,'r-',t5psplit1sthalf(1,:)',meanI1sthalf(1,:)'/startI,'k+','predobs'); hold on;
x1 = t5psplit1sthalf(1,:)';
y1 = meanI1sthalf(1,:)'/startI;
set(fitplot1st,'linewidth',2,'markersize',2);
legend off;

[fitgrph2nd,fit2details] = fit(t5psplit2ndhalf(1,:)',meanI2ndhalf(1,:)'/startI,'exp2')
fitplot2nd = plot(fitgrph2nd,'r-',t5psplit2ndhalf(1,:),meanI2ndhalf(1,:)'/startI,'k+','predobs'); hold on;
x2 = t5psplit2ndhalf(1,:);
y2 = meanI2ndhalf(1,:)'/startI;
set(fitplot2nd,'linewidth',2,'markersize',2);
legend off;

title('Average normalized photobleaching curve and fitting of period 2 to 6')
xlabel('Normalized Time (s)')
ylabel('Normalized I"')
set(gca,'xlim',[0 period],'xtick',[0:period/2:period],'xgrid','on')

%%
%Spectrum change v.s. period
figure(4)
clr = 1;

%Pixel check
if SPE.dim_i < 1024,
  disp('-Invalid data. ROI must be 1024 pixels in horizontal direction.')
  return
end
%Pixel to wavelength data
wav = (1:1024)* 0.50485 + 247.53; %Grating dispersion = 0.50485; Offset = 247.53; http://www.roperscientific.de/gratingcalc.html

%max I / sum laser at max I
Imaxwl = double(IIbase(:,locsls(1,1)+1)) / norm_laserbase(1,locsls(1,1)+1); %Wavelength spectrum at max I / normalized laser at that point
plot(wav,Imaxwl / max(Imaxwl),'color','k'); hold on; %Normalize the max to 1

meanIIline = 1;
for i = locsls(1,1)+2*round(period/(2*nmlfactor*tintvl)):round(period/(2*nmlfactor*tintvl)):locsls(1,1)+(2+9)*round(period/(2*nmlfactor*tintvl)); %Cut framenumber each period
    for j = 1:1:round(period/(2*nmlfactor*tintvl)); %Counting each period
        IIbasehalfp(j,:) = double(IIbase(:,i+j)) / norm_laserbase(1,i+j); %Collect IIbase data and / normalized laser at each point
    end
    
    %Mean of IIbase
    meanIIbasehalfp(meanIIline,:) = mean(IIbasehalfp(:,:)); %Mean of those data
    meanIIline = meanIIline +1;
end

for i = 1:1:10;
    if mod(i,2);
        meanIIbasehalfpplot = plot(wav,meanIIbasehalfp(i,:)/ max(Imaxwl),'color',cc(clr,:),'linewidth',1.5); hold on; %Normalize the Imaxwl to 1
        set(meanIIbasehalfpplot,'LineStyle','-')
    else
        meanIIbasehalfpplot = plot(wav,meanIIbasehalfp(i,:)/ max(Imaxwl),'color',cc(clr,:),'linewidth',1.5); hold on;
        set(meanIIbasehalfpplot,'LineStyle',':')
        clr = clr +1;
    end 
end
set(gca,'xlim',[500 750])
title('Normalized average fluorescence intensity of period 2 to 6')
xlabel('Wavelength (nm)')
ylabel('Normalized average intensity / Laser intensity')
%%
cd(cdpwd)