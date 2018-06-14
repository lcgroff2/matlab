% histogram plot of photon/frame and FWHM

%%
pxlsize = 65.85; %65.85nm/pixel
trkuncerlim = 0.4; % Tracking uncertainty limit in nm
avfwhmlim = 350; %Average FWHM limit in nm
cols = lines(10); 
phminlim = 4; % Photon limit in log10 unit
phmaxlim = 7; % Photon limit in log10 unit
ranpb_s=5:32; % Photobleaching range in sec

%% OPEN FILE

cdpwd = pwd; addpath(pwd);
disp('Choose a folder.')
pathnam = uigetdir('.'); disp(pathnam); cd(pathnam)

fstruct = dir('*.mat');

if length(fstruct)<1,
  disp('-No tracking file (.mat).');
  return;
end

if exist('histph','var')
    clear histph
end

if exist('colcnt','var')
    colcnt = colcnt +1;
else
    colcnt =1;
end

fprintf(1,'Frame number and "average-per-frame" values\n');
for cnt=1:length(fstruct),
    fnam=fstruct(cnt).name;
    load(fnam);
    disp(fnam);
    [numframes,x,ntrack] = size(tracklist{1,1});
    for idx=1:length(tracklist)
        xx=tracklist{idx}(:,1)*pxlsize;
        yy=tracklist{idx}(:,2)*pxlsize;
        fprintf(1,' Track #%i:',idx);
        tmpsz=size(tracklist{idx});
        fprintf(1,' %iframes,',tmpsz(1));
        % # of photon per frame
        avI = round(mean(tracklist{idx}(:,4)));
        stdI = round(std(tracklist{idx}(:,4)))/avI*100; %std in %
        sCMOS_quanteff = 0.57; % This is the quantum efficiency of the sCMOS camera
        collect_eff = 0.03; % This is the collecting efficiency of the microscope 
        conv_gain = 0.6; % Be sure to check ADC setup & camera datasheet for conversion gain
        avph_frame = avI*conv_gain/sCMOS_quanteff/collect_eff;
        % av FWHM
        avfwhm = 2.355*sqrt(mean(tracklist{idx}(:,3)))*pxlsize;
        stdfwhm = 2.355*sqrt(std(tracklist{idx}(:,3)))*pxlsize/avfwhm*100; %std i
        fprintf(1,' I=%4ict, Photon=%1.2e(%1.0f%s), FWHM=%5.1f(%1.0f%s)nm,',avI,avph_frame,stdI,'%',avfwhm,stdfwhm,'%');
        dx=xx(2:end)-xx(1:(end-1));
        dy=yy(2:end)-yy(1:(end-1));
        r_rms=sqrt(mean(dx.^2+dy.^2));
        fprintf(1,' RMSD=%2.1fnm,',r_rms);
        % Tracking uncertainty
        a = 65.85;
        b = 101;
        sigma = avfwhm/(2*sqrt(2*log(2)));
        trkuncer = sqrt(sigma^2/avph_frame+a^2/12/avph_frame+8*pi*sigma^4*b^2/a^2/avph_frame^2);
        fprintf(1,' d=%1.1fnm\n',trkuncer);
        
        ranpb=round(min(ranpb_s)/rnddwelltime_s+1):round(max(ranpb_s)/rnddwelltime_s);
        ranpb_s=ranpb*rnddwelltime_s;
        xpb=msd(xx',ranpb); ypb=msd(yy',ranpb);
        [fitXpb,fitXpbcoef] = fit(ranpb_s',xpb','a*exp(x/(2*T))','startpoint',[1,1]);
        [fitYpb,fitYpbcoef] = fit(ranpb_s',ypb','a*exp(x/(2*T))','startpoint',[1,1]);
        warning('off')
        fitXv = struct(fitXpb); fitXcoefv = struct(fitXpbcoef);
        fitYv = struct(fitYpb); fitYcoefv = struct(fitYpbcoef);
        warning('on')
        ax = fitXv.coeffValues{1,2}; Tx = fitXv.coeffValues{1,1}; R2x_pb = fitXcoefv.rsquare;
        ay = fitYv.coeffValues{1,2}; Ty = fitYv.coeffValues{1,1}; R2y_pb = fitYcoefv.rsquare;
        
        if exist('histph','var')
            histph(end+1,:)=[avI stdI avph_frame avfwhm stdfwhm r_rms trkuncer Tx Ty];
        elseif ~exist('histph','var')
            histph=[avI stdI avph_frame avfwhm stdfwhm r_rms trkuncer Tx Ty];
        end
    end
end
histphtmp=zeros(1,size(histph,2));
for i =1:size(histph,1)
    if histph(i,4)<avfwhmlim
        histphtmp(end+1,:)=histph(i,:);
    end
end
histph = histphtmp(2:end,:);
clear histphtmp

fprintf(1,'"Average-per-frame-per-particle" values\n');
fprintf(1,' %i particles: av FWHM = %3.1f nm with std = %2.1f nm\n',size(histph,1),mean(histph(:,4)),std(histph(:,4)));

for i =1:size(histph,1)
    if histph(i,7) < trkuncerlim
        if exist('histph_l','var')
            histph_l(end+1,:)=histph(i,:);
        elseif ~exist('histph_l','var')
            histph_l=histph(i,:);
        end
    else
        if exist('histph_h','var')
            histph_h(end+1,:)=histph(i,:);
        elseif ~exist('histph_h','var')
            histph_h=histph(i,:);
        end
    end
end
fprintf(1,' %i particles with d<%1.1fnm: av FWHM = %3.1f nm with std = %2.1f nm\n',size(histph_l,1),trkuncerlim,mean(histph_l(:,4)),std(histph_l(:,4)));

figure(1)
ax(1)=subplot(4,1,1);
xvalues=phminlim:(phmaxlim-phminlim)/100:phmaxlim;
[n1, xout1] = hist(log10(histph(:,3)),xvalues);
bar(xout1,n1,'FaceColor', 'none', 'edgecolor', cols(colcnt,:)); grid; hold on;
set(gca, 'xlim', [phminlim phmaxlim])
grid on
ylabel('Occurrence')
title(sprintf('%iptc FWHM=%3.1fnm(%2.1f). %iptc d<%1.1fnm FWHM=%3.1fnm(%2.1f).',size(histph,1),mean(histph(:,4)),std(histph(:,4)),size(histph_l,1),trkuncerlim,mean(histph_l(:,4)),std(histph_l(:,4))))

ax(2)=subplot(4,1,2);
if exist('histph_l','var')
    scatter(log10(histph_l(:,3)),histph_l(:,4),'markerfacecolor', [0.7431 0.8056 0.8056], 'markeredgecolor', cols(colcnt,:)); hold on;
end
if exist('histph_h','var')
    scatter(log10(histph_h(:,3)),histph_h(:,4),'markerfacecolor', 'none', 'markeredgecolor', cols(colcnt,:)); hold on;
end
set(gca, 'xlim', [phminlim phmaxlim],'ylim',[230 350])
box on; grid on;
ylabel('Average FWHM per frame')

ax(3)=subplot(4,1,3);
scatter(log10(histph(:,3)),histph(:,2), 'markeredgecolor', cols(colcnt,:)); hold on;
set(gca, 'xlim', [phminlim phmaxlim], 'ylim', [0 50])
box on; grid on;
ylabel('std of photon per frame')

ax(4)=subplot(4,1,4);
for i =1:size(histph,1)
    if histph(i,8)>1e4
        histph(i,8)=1e4;
    end
    if histph(i,9)>1e4
        histph(i,9)=1e4;
    end
end
scatter(log10(histph(:,3)),(histph(:,8)), 'markeredgecolor', cols(colcnt,:)); hold on;
scatter(log10(histph(:,3)),(histph(:,9)), 'markeredgecolor', cols(colcnt,:)); hold on;
set(gca, 'xlim', [phminlim phmaxlim])
box on; grid on;
ylabel('log(T)')
xlabel('Average photon per frame (log)')

warning('off')
linkaxes(ax,'x');
warning('on')

% save data to a MATfile
fnam2=strtrim(fnam); % remove any spaces
%fnam2=lower(fnam2); % switch to lower case
% find location of '.spe'
fnam2 = strrep(fnam2, '.mat', '');
clear histph_l histph_h
%save('HistPhoton','histph');