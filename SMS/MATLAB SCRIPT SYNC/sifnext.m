
% sifnext:

% Do point kinetics on a set of SIF images and calculate the deathnumber
% for one picked nanoparticle. This version uses previously loaded image
% data, to save time. Also stores key results in "pbkin", a structure.

cc = [0,0,0.75;
    1,0,0;
    0,0.5,0;
    0.75,0.25,0.75;
    0.55,0.27,0.07];

figure(2);
hold on
% imagesc(sumimg'), shading flat % uncomment if you want to sum all frames for image
if curparticle ==1
    imagesc(datacrop.image(:,:,1)'), shading flat; % comment if you want to sum all frames
    set(gca,'dataaspectratio',[1 1 1]) % square pixels
    colorbar
    set(gca,'position',[0.0 0.05 0.9 0.9]); % comment if this makes your figure clip
end
% cax=caxis;
% cax(2)=0.9*cax(2);
% cax(1)=1;
% caxis(cax); % slightly lower one of the caxis values to try to make faint particles visible
 
% note: if image flipped using ', might need to swap xidx, yidx
[xidx,yidx]=ginput(1); %gets N points from the current axes and returns the X- and Y-coordinates in length N vectors X and Y
xidx=round(xidx);
yidx=round(yidx);

% if there is a max one pixel away, change xidx, or yidx
if sum(datacrop.image(xidx-1,yidx,:))>sum(datacrop.image(xidx,yidx,:))
    xidx=xidx-1;
end
if sum(datacrop.image(xidx+1,yidx,:))>sum(datacrop.image(xidx,yidx,:))
    xidx=xidx+1;
end
if sum(datacrop.image(xidx,yidx-1,:))>sum(datacrop.image(xidx,yidx,:))
    yidx=yidx-1;
end
if sum(datacrop.image(xidx,yidx+1,:))>sum(datacrop.image(xidx,yidx,:))
    yidx=yidx+1;
end

if sum(datacrop.image(xidx+1,yidx+1,:))>sum(datacrop.image(xidx,yidx,:))
    xidx=xidx+1;
    yidx=yidx+1;
end
if sum(datacrop.image(xidx+1,yidx-1,:))>sum(datacrop.image(xidx,yidx,:))
    xidx=xidx+1;
    yidx=yidx-1;
end
if sum(datacrop.image(xidx-1,yidx+1,:))>sum(datacrop.image(xidx,yidx,:))
    xidx=xidx-1;
    yidx=yidx+1;
end
if sum(datacrop.image(xidx-1,yidx-1,:))>sum(datacrop.image(xidx,yidx,:))
    xidx=xidx-1;
    yidx=yidx-1;
end

fprintf(1,' First frame pixel intensity: %i\n',datacrop.image(xidx,yidx,1));

pbkin.xpos{curparticle}=xidx;
pbkin.ypos{curparticle}=yidx;
n=5; %% The number of pixels away from the peak. This determines the total number of pixels you want to sum over 
% show markers
for idx = 1:length(pbkin.xpos);
    xpos=pbkin.xpos{idx};
    ypos=pbkin.ypos{idx};
    plot(xpos,ypos,'go')
    rectangle('position',[xpos-n,ypos-n,2*n,2*n],'edgecolor','r')
    txth=text(xpos+n+2,ypos,num2str(idx));
    set(txth,'color',[0 1 0]); % color green
end

% load the rest of file, if this is the first time through
% if curparticle==1,
%     matnam=fnam;
%     matnam(end-2:end)='mat';
%     if exist(matnam,'file') % see if matfile exists
%         disp('mat-file of the same name found, loading...');
%         load(matnam) % load matfile if exists, since this is faster
%     else
%         h=msgbox('loading file, please wait...');
%         %data=sifread([pathnam fnam]); % otherwise just read sif file
%         delete(h);
%     end
% end
    
pkin=zeros([1 data.nframes]);
background=zeros([1 data.nframes]);

%for idx=1:numframes,
%  background(idx)=(data.image(xidx-n,yidx-n,idx)+data.image(xidx-n,yidx+n,idx)+data.image(xidx+n,yidx-n,idx)+data.image(xidx+n,yidx+n,idx))/4;
%  pkin(idx)=sum(sum(data.image((xidx-n):(xidx+n),(yidx-n):(yidx+n),idx)-background(idx)));
%end

% trying to replace for loop with vectorized calculation (need to check!!)
%tic;
background=(datacrop.image(xidx-n,yidx-n,:)+datacrop.image(xidx-n,yidx+n,:)+datacrop.image(xidx+n,yidx-n,:)+datacrop.image(xidx+n,yidx+n,:))/4;
pkin=double(sum(sum(datacrop.image((xidx-n):(xidx+n),(yidx-n):(yidx+n),:))))-(2*n+1)^2*double(background);
pkin=squeeze(pkin); % removes singleton dimensions
%toc

sCMOS_quanteff = 0.57; % This is the quantum efficiency of the sCMOS camera
collect_eff = 0.03; % This is the collecting efficiency of the microscope 
conv_gain = 0.6; % Be sure to check ADC setup & camera datasheet for conversion gain
pkin_frame=pkin*conv_gain/sCMOS_quanteff/collect_eff;
deathnumber=sum(pkin_frame);
photonsperframe = deathnumber/datacrop.nframes;

figure(3);
hold all
plot(pkin,'color',cc(curparticle,:));

legendText{curparticle} = sprintf('Particle %i',curparticle);
legend(legendText{:});
tstr1=sprintf('Intensity track of file: %s\n',fnam);
%tstr1=sprintf('Intensity track of file: %s\n nanoparticle #%i, at x= %i and y= %i\n with total photons = %1.5g',fnam,curparticle,xidx,yidx,deathnumber);
title(tstr1)
xlabel('Frame number');
ylabel('Intensity');
set(gca,'xlim',[0 endf]);
%figure (4);
%plot(pkin/max(pkin));
%tstr=sprintf('Intensity track of file: %s\n the picked nanoparticle locates at x= %3.1f and y= %3.1f',fnam,xidx,yidx);
%title(tstr);

% quick estimate of 1/e time by seeing where trace crosses 1/e
%[minval,minidx]=min(abs(pkin-pkin(1)*0.632));
%tau_est = minidx*dwelltime_s;
% estimate tau using integral (sum) method
Iint=sum(pkin);
tau_est=sum(pkin(:)'.*datacrop.xdat)/Iint;

% fraction of total decay
frac_decay = 1-pkin(end)/pkin(1);

fprintf(' Estimated total photons emitted before death (deathnumber) = %1.5g\n',deathnumber);
fprintf(' Estimated photons per frame) = %i\n',photonsperframe);
fprintf(1, ' Particle position(x,y): %i, %i \n', xidx, yidx);
fprintf(1,' Decay roughly %5.1f percent complete \n',100*frac_decay);
fprintf(1,' Rough estimate of 1/e time: %1.4g s\n',tau_est);

% save data into a structure/list, to make it easier to summarize results
% or analyze multiple kinetics traces
pbkin.kin{curparticle}=pkin;
pbkin.tau{curparticle}=tau_est;
pbkin.deathnum{curparticle}=deathnumber;
pbkin.frac_decay{curparticle} = frac_decay;
fprintf(1,' This kinetics trace saved as pk.kin{1,%i}\n',curparticle);

