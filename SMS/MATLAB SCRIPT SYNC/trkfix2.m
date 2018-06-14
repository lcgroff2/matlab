
% TRKFIX -- version 7 of SPE tracking script, with gaussian fit

% this an "automatic" version is for running "in the background", analyzing
% several trajectories automatically without needing to type in the
% names.

% do not run this script directly -- instead use "trkfixrun"

% Note: this is for "fixed" particles

% TODO: 
% -need better "corrplot" that asks for how many frames, first
% makes a summary of "all" important parameters for "all"
% trajectories, and then asks you to "pick" which one to plot.

%DONE:
% -save as .mat file with same name but .mat extension instead of
% .spe
% -spetrkrun will define variables that define how threshold is set
% -use the "raw" particle position data, without running the
% tracking function.
% -only save important variables: tracklist, a summed image, etc.
% -better background subtraction (uses the min of the median)
% -better automatic threshold calculation.
% -if possible, only use brightest 2-3 particles for vibration
% correction.
% -find peaks based on summed image, do pruning,
% then use the *same* initial peak positions for fitting *all* of the frames.
% -use my "fpass" instead of "bpass"

% params to set in trkfixrun:
%   bp1, bp2: first & 2nd bandpass parameters)
%   pksep: min separation between peaks

% Note: this version does not use Weeks' "track" function, but
% instead keeps the same order of peak locations given by the
% "peaks" function.  To revert to Weeks' code, uncomment the
% indicated line near the end of this script.

fprintf(1,'\n%i/%i\n',cnt,length(fstruct));
fprintf(1,'FILE\n');
fprintf(1,'%s\n',fnam);

if ~isempty(regexp(fnam, '.SPE', 'once')) || ~isempty(regexp(fnam, '.spe', 'once'))
    data = loadspe(fnam);
elseif ~isempty(regexp(fnam, '.sif', 'once'))
    data = sifread(fnam);
end

if ~isempty(regexp(fnam, '.sif', 'once'))
    %% CROP
    fprintf(1,'CROPPED DATA\n');
    data.rnddwelltime_s = round(data.dwelltime_s*1000)/1000;
    Ibright = data.image(data.xmaxI,data.ymaxI,:);
    Ibright = Ibright(1,:);
    diffIbright = double([diff(Ibright),0]);
    frame = 1:data.nframes;
    figure(2)
    [ax,Ibr] = plotyy(frame,Ibright,frame,diffIbright); 
    set(ax,'xlim',[0 data.nframes]);
    set(get(ax(1),'ylabel'),'string','I');
    set(get(ax(2),'ylabel'),'string','Diff(I)');
    xlabel('Framenumber')
    title('Monitoring brightest pixel in order to crop frames')
    [mx,start] = max(diffIbright);
    if mean(data.image(data.xmaxI,data.ymaxI,1:start)) > 125;
        i = 2;
        while data.image(data.xmaxI,data.ymaxI,i) < 1.5*mean(data.image(data.xmaxI,data.ymaxI,1:i-1)) && i < 2/data.rnddwelltime_s % Find a threshold within first 2 seconds
            i = i + 1;
        end
        start = i - 1;
    end
    startf = start +1;
    endf = start + round(cropsec/data.rnddwelltime_s); % Crop ends at 'cropsec'
    datacrop.nframes = endf - startf +1;
    datacrop.xdat = ((1:datacrop.nframes)-1)*data.dwelltime_s;
    datacrop.image = data.image(:,:,startf:endf);
    fprintf(1,' Cropped size: %ix%ix%i\n',data.width,data.height,datacrop.nframes);
    fprintf(1,' Frame number was cropped from #%i to #%i\n',startf,endf);
    if startf - 3 > 0
        fprintf(1,' I at #Frame number of the brightest particle (*starting crop data):\n');
        I_3 = data.image(data.xmaxI,data.ymaxI,startf-3);
        I_2 = data.image(data.xmaxI,data.ymaxI,startf-2);
        I_1 = data.image(data.xmaxI,data.ymaxI,startf-1);
        I0 = data.image(data.xmaxI,data.ymaxI,startf);
        I1 = data.image(data.xmaxI,data.ymaxI,startf+1);
        I2 = data.image(data.xmaxI,data.ymaxI,startf+2);
        I3 = data.image(data.xmaxI,data.ymaxI,startf+3);
        fprintf(1,' I#%i=%i I#%i=%i I#%i=%i *I#%i=%i I#%i=%i I#%i=%i I#%i=%i\n',startf-3,I_3,startf-2,I_2,startf-1,I_1,startf,I0,startf+1,I1,startf+2,I2,startf+3,I3);
    end
    datareal = data;
    data = datacrop;
end

%%
[dim_i,dim_j,numframes] = size(data.image);

% array for storing all the centers
centerlist=[];

imghandle=[];
cross=[];

firstimg=data.image(:,:,1);

sumimg=sum(data.image,3); % sum up all frames (along third dimension,
                    % which is frame number)
% Make fourier filter function
H=mkffilt(sumimg,bp1,bp2);
% use Jason's fourier filter function 
bpimg=fpass(sumimg,H);

% find particle initial approximate positions only once, but use centerg to
% find the subpixel position in each frame.
imgmax=max(max(bpimg));
imgsz=size(bpimg);
imgmin=min(min(bpimg((bp2+1):(imgsz(1)-bp2-1),(bp2+1):(imgsz(2)-bp2-1))));
thresh=imgmin+(imgmax-imgmin)*threshfrac;
peaks=pkfnd(bpimg,thresh,minsep);
peaksz=size(peaks);
centers=zeros([peaksz(1) 4]);

% kludge to decrease thresh until number of particles is ~20
while peaksz(1)>maxnumparticles,
  thresh=thresh+0.5*thresh;
  peaks=pkfnd(bpimg,thresh,minsep);
  peaksz=size(peaks);
end



% do gaussian fitting to summed data to refine peak positions
for idx=1:peaksz(1), % loop over multiple peaks
    centers(idx,:)=cntr3dg(sumimg,peaks(idx,:),gpixels); % only do cntr3dg on raw img
    peaks(idx,:)=centers(idx,1:2);
end

% check for "duplicates" -- peaks that, after refinement, are at
% same position.
for idx=1:(peaksz(1)*3), % *3 is experimental kludge, shouldn't
                         % hurt anything
  for idx2=1:peaksz(1),
    xd=(peaks(idx2,1)-peaks(:,1)).^2;
    yd=(peaks(idx2,2)-peaks(:,2)).^2;
    r=sqrt(xd+yd);
    dupidx=find(r<minsep); % a dupe if separation is too small
    tmpidx=find(dupidx ~= idx2);
    if ~isempty(tmpidx),
      dupidx=dupidx(tmpidx);
      % remove the last duplicate
      dupidx=dupidx(end);
      if dupidx==1, % shouldn't happen, but include just in case
	peaks=peaks(2:peaksz(1),:);
      elseif dupidx==peaksz(1),
	peaks=peaks(1:(peaksz(1)-1),:);
      else
	peaks=peaks([1:(dupidx-1) (dupidx+1):peaksz(1)],:);
      end
      % recalculate peaksz
      peaksz=size(peaks);
      break
    end
  end
end

unsorttrk=zeros([numframes peaksz(1) 4]);
centers=zeros([peaksz(1) 4]);

%%
cc = lines(25);
for idx=1:numframes,
    img=double(data.image(:,:,idx));
    % analysis method based on tutorial at
    % http://physics.georgetown.edu/matlab/tutorial.html
    warning off
    if ~isempty(peaks),
        for idx2=1:peaksz(1), % loop over multiple peaks
            centers(idx2,:)=cntr3dg(img,peaks(idx2,:),gpixels); % only do centerg on raw img
	        %centers(idx2,:)=centerg(img,peaks(idx2,:),gpixels);
        end
    else
        centers=[];
    end
    warning on
    centersz=size(centers);
    % store "unsorted" positions
    if ~isempty(centers),
      unsorttrk(idx,:,:)=centers;
    end
    % do some plotting if showimages==1
    if showimages==1 && idx~=numframes
        figure(1)
        if isempty(imghandle),
            imagesc(img);
            %colorbar
            %colormap('jet');
            imghandle=get(gca,'children');
            set(gca,'YDir','normal')
        else
            set(imghandle,'cdata',img);
            if ~isempty(cross);
                delete(cross);
                cross=[];
                if ~isempty(txth);
                    delete(txth);
                end
                txth=[];
            end
            minI=min(min(img));
            maxI=max(max(img));
        end
        if ~isempty(centers),
            for idx2=1:centersz(1),
                xx=centers(idx2,1);
                yy=centers(idx2,2);
                cross(idx2)=line([xx-4 xx+4 xx xx xx],[yy yy yy yy-4 yy+4]);
                set(cross(idx2),'color',[0 1 0]);
                txth(idx2)=text(xx+2,yy-4,num2str(idx2));
                set(txth(idx2),'color',[0 1 0]); % color green
                %legendText{idx2} = sprintf('%i',idx2);
            end
        end
	  drawnow; % make sure the image updates
      %legend(legendText{:});
      %pause(imagedelay);
    elseif (showimages==0 && idx==1) || (showimages==1 && idx==numframes)
        img=double(data.image(:,:,1));
        figure(1)
        if isempty(imghandle),
            imagesc(img);
            colorbar
            colormap('jet');
            imghandle=get(gca,'children');
            set(gca,'YDir','normal')
        else
            set(imghandle,'cdata',img);
            colorbar
            colormap('jet');
            if ~isempty(cross);
                delete(cross);
                cross=[];
                if ~isempty(txth);
                    delete(txth);
                end
                txth=[];
            end
            minI=min(min(img));
            maxI=max(max(img));
        end
        if ~isempty(centers),
            for idx2=1:centersz(1),
                xx=centers(idx2,1);
                yy=centers(idx2,2);
                cross(idx2)=line([xx-4 xx+4 xx xx xx],[yy yy yy yy-4 yy+4]);
                set(cross(idx2),'color',[0 1 0]);
                txth(idx2)=text(xx+2,yy-4,num2str(idx2));
                set(txth(idx2),'color',[0 1 0]); % color green
            end
        end
        drawnow; % make sure the image updates
    end
	
    if ~isempty(centers),
        tmpcenterlist=centers;
        % use the time index as the third column
        tmpcenterlist(:,5)=idx;
        % now concatenate onto our centerlist
        if length(centerlist)<1, % special case for first center
            centerlist=tmpcenterlist;
        else
            centerlist=cat(1,centerlist,tmpcenterlist);
        end
    end
end

%%
% Now take the positions and times and feed to the track.m function
% columns 1:5 in centerlist are X, Y, intensity, square spotsize,
% and Time (frame number), respectively.
% but positionlist muxt be in the order: X, Y, other1, other2, time
% so I rearrange here:
if length(centerlist)<1,
  fprintf(1,'Bad data set: %s\n',fnam);
  cd(cdpwd)
  return
end
positionlist=centerlist(:,[1 2 3 4 5]); 
trackparam.mem=1;
trackparam.good=200;
trackparam.dim=2;
trackparam.quiet=0;
%tracks=track(positionlist,5,trackparam);

% convert array to a list of arrays called "tracklist":
%numtracks=tracks(end,end);
%tracklist={}; % create an empty list

%for idx=1:numtracks,
%    curtrackidx=find(tracks(:,6)==idx);
%    tmplist=tracks(curtrackidx,1:5);
%    tracklist{idx}=tmplist;
%end

unsorttrksav=unsorttrk; % for debugging only

% remove duplicates from unsorted tracks, looking at first point
for idx=1:(peaksz(1)*3), % *3 is experimental kludge, shouldn't
                         % hurt anything
  for idx2=1:peaksz(1),
    xd=(unsorttrk(1,idx2,1)-unsorttrk(1,:,1)).^2;
    yd=(unsorttrk(1,idx2,2)-unsorttrk(1,:,2)).^2;
    r=sqrt(xd+yd);
    dupidx=find(r<minsep); % a dupe if separation is too small
    tmpidx=find(dupidx ~= idx2);
    if ~isempty(tmpidx),
      dupidx=dupidx(tmpidx);
      % remove the last duplicate
      dupidx=dupidx(end);
      if dupidx==1, % shouldn't happen, but include just in case
	unsorttrk=unsorttrk(:,2:peaksz(1),:);
      elseif dupidx==peaksz(1),
	unsorttrk=unsorttrk(:,1:(peaksz(1)-1),:);
      else
	unsorttrk=unsorttrk(:,[1:(dupidx-1) (dupidx+1):peaksz(1)],:);
      end
      % recalculate peaksz
      peaksz=size(unsorttrk);
      peaksz=[peaksz(2) 2]; % kludge so I could cut & paste from above
      break
    end
  end
end

clear unsrttracklist

% record unsorted tracks
for idx=1:peaksz(1),
  unsrttracklist{idx}=unsorttrk(:,idx,:);
end



% now the list of trajectories is in tracklist{particlenumber}
% to get all the "x" trajectories for the first particle trajectory,
% specify: tracklist{1}(:,1)
% for the "y" trajectories for the first particle, 
% specify: tracklist{1}(:,2)
% for brightness of first particle, specify:
% tracklist{1}(:,3)
% for the variance (width squared), is "4"

% NOTE: be careful not to confuse curly braces "{}" with parentheses "()"



% clean up workspace a bit
clear txth cross idx idx2 tmpcenterlist centersz trackparam %centers 
clear curtrackidx curxtrack curytrack numtracks thresh xx yy curcolor colorstr
clear tmplist imghandle fitlen dx dy peaks peaksz centerlist dim_i dim_j positionlist bpimg r_rms

% save data to a MATfile
fnam2=strtrim(fnam); % remove any spaces
%fnam2=lower(fnam2); % switch to lower case
% find location of '.spe'
fnam2 = strrep(fnam2, '.SPE', ''); fnam2 = strrep(fnam2, '.sif', ''); fnam2 = strrep(fnam2, '.', '_');

% comment next line to use the Weeks track routine instead.
tracklist=unsrttracklist;

fprintf(1,'Got %i tracks\n',length(tracklist));
rnddwelltime_s = datareal.rnddwelltime_s;
save(fnam2,'tracklist','firstimg','fnam','rnddwelltime_s');
saveas(figure(1),[fnam2 '-IMG']);
%saveas(figure(2),[fnam2 '-Ibrightest']);

