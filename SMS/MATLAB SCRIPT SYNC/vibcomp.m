
% "vibcomp.m" -- vibration compensation script.  Only use for
% "stationary" particles.  Not for brownian motion in solution.

% Gives corrected tracklist in list variable "corrtracklist".

% set some threshold values for judging the trajectories.  You may
% adjust these if you feel some valid trajectories have been
% rejected.

maxfwhm = 300;
%maxrmsd = 16; %=32.9 nm/frame
%minint = 1000;
% maxfwhm=5 pixel; =329.25nm
% maxrmsd=0.07 pixel/frame; =4.61 nm/frame
% minint=1000;

% Uncomment next line if you want to only analyze the first 500 points.
% Or uncomment and set to "1000" to analyze only the first 1000
% points, etc.
%truncate=1000;

tmptracklist=tracklist;

% End of user-configurable parameters.  Hopefully no adjustment is
% required beyond this line. 

% Do truncation if "truncate" variable is set.
if exist('truncate')==1,
  fprintf(1,'Truncating to first %i points.\n',truncate);
  for idx=1:length(tmptracklist),
    tmptrk=tmptracklist{idx};
    tmpsz=size(tmptrk);
    if tmpsz(1)>truncate,
      tmptracklist{idx}=tmptrk(1:truncate,:);
    end
  end
end


% make some variables to "judge" the trajectories
trklen=zeros(size(tmptracklist));
trkint=trklen;
trkfwhm=trklen;
trkrmsd=trklen;

for idx=1:length(tmptracklist)
    xx=tmptracklist{idx}(:,1)*pxlsize;
    yy=tmptracklist{idx}(:,2)*pxlsize;
    fprintf(1,' Track #%i:',idx);
    tmpsz=size(tmptracklist{idx});
    fprintf(1,' frame = %i,',tmpsz(1));
    trklen(idx)=tmpsz(1);
    fprintf(1,' av I = %i, av FWHM = %5.1f nm,',round(mean(tmptracklist{idx}(:,4))),2.355*sqrt(mean(tmptracklist{idx}(:,3)))*pxlsize);
    trkint(idx)=round(mean(tmptracklist{idx}(:,4)));
    trkfwhm(idx)=2.355*sqrt(mean(tmptracklist{idx}(:,3)));
    dx=xx(2:end)-xx(1:(end-1));
    dy=yy(2:end)-yy(1:(end-1));
    r_rms=sqrt(mean(dx.^2+dy.^2));
    trkrmsd(idx)=r_rms;
    fprintf(1,' av RMSD = %2.1f nm/frame\n',r_rms);
end

% now judge which tracks are good enough
rmsdlimit = input('RMSD maximum limit = ');
Ilimit = input('I minimum limit = ');
goodtracks=[];
maxtrklen=max(trklen);

for idx=1:length(tracklist),
  if (trkint(idx)>=Ilimit) && (trklen(idx)==maxtrklen) && ...
	(trkfwhm(idx)<=maxfwhm) && (trkrmsd(idx)<=rmsdlimit),
    % append list of good tracks
    goodtracks=[goodtracks idx];
  end
end

if length(goodtracks)<2,
  fprintf(1,'Insufficient good tracks found. Check criteria.\n');
  return;
else
  fprintf(1,'Number of good tracks found is %i. Tracks: ', length(goodtracks));
  fprintf(1,'%i ',goodtracks);
  fprintf(1,'\n');
end

% select only "good" tracks
clear corrtracklist
for idx=1:length(goodtracks),
  corrtracklist{idx}=tmptracklist{goodtracks(idx)};
end

% calculate and apply vibration correction to corrtracklist

% some temporary variables to store mean positions
meanx=0;
meany=0;

for idx=1:length(corrtracklist),
  meanx=meanx+corrtracklist{idx}(:,1);
  meany=meany+corrtracklist{idx}(:,2);
end

% divide by number of tracks
meanx=meanx/length(corrtracklist); 
meany=meany/length(corrtracklist);

% recenter around zero
meanx=meanx-mean(meanx);
meany=meany-mean(meany);

% subtract meanx, meany from each of the trajectories
for idx=1:length(corrtracklist),
  corrtracklist{idx}(:,1)=corrtracklist{idx}(:,1)-meanx;
  corrtracklist{idx}(:,2)=corrtracklist{idx}(:,2)-meany;
end

fprintf(1,'Corrected RMSD for the good tracks\n');
% recalculate rmsd for comparison to uncorrected trajectories
for idx=1:length(corrtracklist)
    xx=corrtracklist{idx}(:,1)*pxlsize;
    yy=corrtracklist{idx}(:,2)*pxlsize;
    fprintf(1,' Track #%i,',goodtracks(idx));
    dx=xx(2:end)-xx(1:(end-1));
    dy=yy(2:end)-yy(1:(end-1));
    r_rms=sqrt(mean(dx.^2+dy.^2));
    fprintf(1,' av RMSD = %2.1f nm/frame\n',r_rms);
end
    
% clear some temporary variables
%clear trklen trkint trkfwhm truncate tmpsz tmptrk maxfwhm maxrmsd ...
%    goodtracks
%clear tmptracklist r_rms meanx meany xx yy dx dy;
