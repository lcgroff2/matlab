
% sifkin.m
% Do point kinetics on a set of SIF images and calculate the deathnumber
% for picked nanoparticles. Uses sifnext.m to get the next particle, and sifreport.m
% to generate a brief report. Draws a box indicating the summed area for each particle,
% and numbers each particle.

% Stores results in "pk", a list of structures with the following elements, where nfile
% is the file number (number increments each time a file is loaded), and nparticle is
% the particle number.
%pk{nfile}.name = file name
%pk{nfile}.xpos{nparticle} = x pos of particle
%pk{nfile}.ypos{nparticle} = y pos of particle
%pk{nfile}.xdat = time axis (same for every kinetics trace for this file)
%pk{nfile}.ydat{nparticle} = kinetics trace data for particle
%pk{nfile}.deathnum{nparticle} = death number of particle
%pk{nfile}.tau{nparticle} = estimated lifetime

% example: to plot kinetics from the fifth particle from the third file analyzed this session:
%>> plot(pk{3}.xdat,pk{3}.ydat{5})
%>> title(['file: ' pk{3}.name ])

% NOTE 1: if you want to be sure to include very faint particles in your analysis,
% then replace 'imagesc' with 'logimagesc', and comment out 'colorbar', since
% 'logimagesc' makes its own log-scaled colorbar

% NOTE 2: If you want to use the "summed" image, then edit the imagesc lines as per comments

% wish list:
% - make biexponential/kww fitting routine, that picks file numbers, etc.
% - make autocorrelation routine
% - show all 3 image types and give option of which to use (?)
% - make "sifshow" script
% - make particle tracking (centroid) script
% - option to reset the color axis, ignoring brightest particles (just for color purposes)

% done:
% - for each file, make pk{1}, pk{2}, etc. with all important results
% - make sifread convert image to uint16
% - try to find nearby peaks, and re-center
% - make sifread1st, to get first frame, then read in background? to make analysis faster
% - pre-load and save as a 16-bit mat file, with same file name, to speed up.

close all
clear

if exist('pbkin','var')
    %disp('warning: will overwrite pbkin variable.')
  clear pbkin % no warning needed: data saved in pk list/struct
end

%% OPEN FILE
if ispc || ismac, %ISPC--True for the PC (Windows) version of MATLAB.ISPC, returns 1 for PC (Windows) versions of MATLAB and 0 otherwise.
  [fnam, pathnam, filterindex] = uigetfile('*.sif', 'Pick an SIF-file');
  %cd(pathnam); 
else
  ls *.sif
  fnam=input('enter SIF file name ','s');
end
disp('FILE')
fprintf(1,'%s\n\n',fnam);

%% 
%h=msgbox('loading file, please wait...');
data = sifread([pathnam fnam]); % Read frames
%delete(h)
%% CROP
fprintf(1,'\nCROPPED DATA\n');
Ibright = data.image(data.xmaxI,data.ymaxI,:);
Ibright = Ibright(1,:);
diffIbright = double(abs([diff(Ibright),0]));
frame = 1:data.nframes;
figure(1)
[ax,Ibr] = plotyy(frame,Ibright,frame,diffIbright); 
set(ax,'xlim',[0 data.nframes]);
set(get(ax(1),'ylabel'),'string','I');
set(get(ax(2),'ylabel'),'string','Diff(I)');
xlabel('Framenumber')
title('Monitoring brightest pixel in order to crop frames')
[mx,start] = max(diffIbright);
if mean(data.image(data.xmaxI,data.ymaxI,1:start)) > 125;
    start = 0;
end
startf = start +1;
endf = start + round(5/data.dwelltime_s); % Crop ends at 160 sec
datacrop.nframes = endf - startf +1;
datacrop.xdat = ((1:datacrop.nframes)-1)*data.dwelltime_s;
datacrop.image = data.image(:,:,startf:endf);
fprintf(1,'Cropped size: %ix%ix%i\n',data.width,data.height,datacrop.nframes);
fprintf(1,'Frame number was cropped from #%i to #%i\n',startf,endf);

%%
%[dim_i,dim_j] = size(data.image(:,:,1));
%numframes=data.nframes;

%%

%sumimg=sum(data.image,3); % use for plotting, in case some particles are "off" at first frame,
                          % or first frame is blank due to shutter timing slightly slow
figure(2);
hold off
% imagesc(sumimg'), shading flat % uncomment if you want to sum all frames for image
imagesc(datacrop.image(:,:,1)'); % comment if you want to sum all frames
set(gca,'dataaspectratio',[1 1 1])
colorbar
set(gca,'position',[0.0 0.05 0.9 0.9]); % comment if this makes your figure clip
drawnow % force the plot to show before continuing

curparticle=1;

%% loop starts here
fprintf(1,'\nDATA AT SELECTED PARTICLES\n');
while 1 % keep looping until "break" reached
  if curparticle==1
    button=questdlg('Image OK, want to analyze?');
  else
    button=questdlg('Analyze another particle?');
  end
  if button(1)=='Y',
    fprintf(1,'Particle %i\n',curparticle);
    sifnext;
    curparticle=curparticle+1;
  elseif button(1)=='N' && curparticle==1
    disp('No particle selected')
    break;
  else
    break; % exit the loopc
  end
end

%% do report summary
if exist('pbkin','var')
  fprintf(1,'\nSUMMARY REPORT\n');
    for idx = 1:length(pbkin.xpos);
      fprintf(1,'Particle %i: Deathnum = %1.4g, Tau_est = %1.4g s, Decay %2.3g pct\n',idx,pbkin.deathnum{idx},pbkin.tau{idx},pbkin.frac_decay{idx}*100);
    end
%   if exist('pk','var')
%       nfile = length(pk)+1;
%   else
%       nfile = 1;
%   end
  % now save results into an "overall" struct, one entry per sif file
  pk.name = fnam;
  pk.date = datestr(now,'mm-dd-yy HH:MM PM');
  pk = pbkin;
  pk.image = datacrop.image(:,:,1); % just save first frame
  savefile=questdlg('Want to save the data and figures?');
  if savefile(1) == 'Y'
     directory = [datestr(now,'mm-dd-yy HH.MM.SS ') fnam];
     mkdir(directory)
     save([directory '/pk'], 'pk');
     saveas(figure(2),[directory '/figure2']);
     saveas(figure(3),[directory '/figure3']);
  end
end

%% now check for duplicates

% if exist('pk','var')
%     for idx=1:(length(pk)-1)
%         if strcmp(pk{idx}.name,pk{nfile}.name)
%             resp=questdlg('Previous analysis exists on this file. Select action','Duplicate','keep both','delete older analysis','delete current analysis');
%             switch resp,
%               case 'keep both'
%                 disp('keeping both');
%               case 'delete older analysis'
%                 pk(idx)=[];
%               case 'delete current analysis'
%                 pk(nfile)=[];
%             end
%         end
%     end
% end

%%clear pbkin