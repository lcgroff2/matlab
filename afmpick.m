
% AFMPICK -- script for plotting/analyzing AFM results/give a histogram and
% a gaussian distribution curve

% To use, first change to the directory/folder containing the AFM file,
% then run AFMPICK.  Then select the AFM file.  It uses LOADAFM to load
% the data.  Then use left mouse to select particles.  The script
% will locate the nearest local maximum, so you don't need to be perfectly
% accurate.  The line scan, going through the local maximum, is shown in
% the right window.  Particle height is calculated as the peak value minus the 
% averages of the two baseline values on either side of peak.  If you click
% on a "bad" particle, you can "right-click" to remove the most recent
% particle from the list of particles.  When you have enough particles for
% your histogram, then click the mouse wheel to end.  The particle heights
% are stored in the variable "plistz".

% This script was developed/tested on an image of roughly 15-50 nm height particles,
% with 10 nm pixels (5 microns scan size, 500 lines). 
% also tested for the big particles (~250 nm). It is OK, however,you may
% change the histgoram with "hwidth" according to the different big
% particle size.
%  If you use a different pixel size, some parameters in the code must be changed -- send
% me a datafile and I can fix it in a couple of minutes.

% -JDM 6/2010

% histogram and pub-style AFM image added -JY 6/2010

if ispc,
  [fnam, pathnam, filterindex] = uigetfile('*.afm', 'Pick an AFM-file');
  cd(pathnam);
else
  ls *.afm
  fnam=input('enter AFM file name ','s');
end

afm=loadafm(fnam);

if isempty(afm),
  return;
end

figure(1)
delete(1) % clear the figure, else figure might have strange props that make this script not work properly
figure(1)
set(1,'units','pixels');
figpos=get(1,'position');
% rect = [left, bottom, width, height]
scrsz=get(0,'screensize');

% set figure to be 1300x700 pixels, if user screen large enough
if scrsz(3)>=1350,
  set(1,'position',[5 45 1300 700]);
else
  set(1,'position',[5 25 1000 600]);
end

Z=afm.zdata;
XY=afm.xy_range;
nl=afm.image_res; % number of lines
fprintf(1,'File name: %s\n',fnam);
fprintf(1,'Description: %s\n',afm.short_desc);
fprintf(1,'Plot size %3.1f microns\n',afm.real_size);
fprintf(1,'Number of lines: %i\n',afm.image_res);
fprintf(1,'Pixel size: %5.2f nm\n',afm.real_size/afm.image_res*1000);

ax1=subplot('position',[.02 .1 .45 .75]);

% Z=Z(:,end:-1:1);
Z=Z(end:-1:1,:);
imagesc(XY,XY,Z)
%pause
axis ij % sets origin in upper left corner
axis square % makes pixels square
shading interp
colormap hot
colorbar

ax2=subplot('position',[.55 .2 .4 .6]);

pixsz=afm.real_size/afm.image_res;
nl=afm.image_res;

fprintf(1,'\n\n  *** Click left button to pick a particle, right button to\n delete last/previous particle, click wheel to finish ***\n\n');

set(1,'currentaxes',ax1);
hold on
heightlist=0;
button=0;
idx=1;
plistx=[];
plisty=[];
plistz=[];
while button ~= 2,
  set(1,'currentaxes',ax1);
  [x,y,button]=ginput(1);
  if button==2, % exit point selection loop
    break;
  end
  if button==3, % remove most recent point
    if txth==[],
      break;
    end
    idx=idx-1;
    delete(txth);
    txth=[];
    plistx=plistx(1:end-1);
    plisty=plisty(1:end-1);
    plistz=plistz(1:end-1);
    continue;
  end
  % convert to pixel units;
  x=round(x/pixsz)+1;
  y=round(y/pixsz)+1;
  if x<3,
    x=3;
  end
  if y<3,
    y=3;
  end
  if x>(nl-3),
    x=nl-3;
  end
  if y>(nl-3),
    y=nl-3;
  end
  % refine peak position
  rn=-3:3;
  [mx,idx1]=max(max(Z(y+rn,x+rn)));
  [mx,idx2]=max(max(Z(y+rn,x+rn)'));
  nx=x+idx1-4;
  ny=y+idx2-4;
  %idx1
  %idx2
  %Z(y+rn,x+rn)
  if Z(ny,nx)>=Z(y,x),
    %disp('refinement worked')
    x=nx;
    y=ny;
  else
    disp('BZZT.  Need to swap x,y')
  end
  hold on
  plot((x-1)*pixsz,(y-1)*pixsz,'go')
  txth=text((x+1.5)*pixsz,(y-2)*pixsz,num2str(idx));
  set(txth,'color',[0 1 0]);
  
  % now plot data
  ll=x-2;
  ul=x+2;
  if ll<2,
    ll=2;
  end
  if ul>(nl-1)
    ul=nl-1;
  end
  % extend the range until local minimum is found on either side
  while (Z(y,ll)>Z(y,ll-1)) & (Z(y,ul)>Z(y,ul+1));
    ll=ll-1;
    ul=ul+1;
    if ll<=1,
      ll=2;
      break;
    end
    if ul>=(nl-1),
      ul=nl-1;
      break;
    end
  end
  % add three more points
  if ll>4,
    ll=ll-3;
  end
  if ul<(nl-4),
    ul=ul+3;
  end
  
  rn=ll:ul;
  set(1,'currentaxes',ax2);
  hold off
  plot(XY(rn),Z(y,rn));
  % fix axes up a bit
  ax=axis;
  ax(1)=min(XY(rn));
  ax(2)=max(XY(rn));
  ax(3)=min(Z(y,rn));
  ax(4)=max(Z(y,rn));
  ax(3)=ax(3)-.03*(ax(4)-ax(3));
  ax(4)=ax(4)+.03*(ax(4)-ax(3));
  axis(ax);
  hold on
  % height of plane on left side=min
  plmin=min(Z(y,ll:(ll+3)));
  % plane height on right side
  prmin=min(Z(y,(ul-3):ul));
  % peak height
  ppeak=Z(y,x);
  abspeak=ppeak-(plmin+prmin)/2;
  % plot line showing peak,
  plot([XY(x) XY(x)],[ppeak (plmin+prmin)/2],'r')
  text(XY(x),(ppeak+plmin)/2,num2str(abspeak,3));
  hold off
  plistx(idx)=XY(x);
  plisty(idx)=XY(y);
  plistz(idx)=abspeak;
  idx=idx+1;
end

hold off

% display results,

fprintf(1,'particle results:\n');
fprintf(1,' part #     height     xpos     ypos\n');
for idx=1:length(plistx),
  fprintf(1,'  %2i         %4.1f     %5.2f    %5.2f\n', idx, plistz(idx), plistx(idx), plisty(idx))
end

fprintf(' particle heights in <plistz> variable\n');
%plot an AFM image for publication
% delete (1);
% figure (2); 
% subplot(1,2,1);
figure(2)
imagesc(XY,XY,Z)
%pause
axis ij % sets origin in upper left corner
axis square % makes pixels square
% axis off
shading interp
colormap hot
colorbar
text(6,1.8,'Nanoparticle Height (nm)','Rotation',-90,'FontSize',12);
line([0.4 0.9],[4.8 4.8],'LineStyle','-','Color','w','LineWidth',4);
text(0.3,4.6,'500 nm','Color','w','FontSize',14);
axis off
figure(1)

% histgauss.m -- make a histogram, with overlaid gaussian
% edit to change the number and width of bins

% assume "h" is the variable that contains height data

h=plistz;
disp('assuming h is the heights variable');

if length(h)<50,
    disp('not enough height values');
end

% start value for histogram
hstart=0;

% width of each bin (1, 2, 2.5, 3, or 5 are reasonable values)
hwidth=2.5;

% end value for histogram
hend=max(h);

edges=hstart:hwidth:hend;

% calculate histogram, assuming "h" is height data
N=histc(h,edges);

% plot histogram.  How do I know this command works?  I typed "help hist",
% which led me to "histc" which seemed better because you could specify
% edges.  "help histc" gives the plotting command.
% subplot(1,2,2);

set(1,'currentaxes',ax2);
bar(edges,N,'histc');
hbar = findobj(gca,'Type','patch');
set(hbar,'FaceColor','b','EdgeColor','k');

% Now calculate the gaussian center and width, but first remove some that are
% outliers.  No "one correct way" to do this, but this seems more or less
% reasonable:
% Usually, it is reasonable to to remove 5% smallest and biggest height
% values and then do the the gaussian distribution analysis

h=sort(h);
n_omit=floor(length(h)*.05);
if n_omit<2,
 n_omit=1;
end

h=h(n_omit:(end-n_omit));


% find all the heights between the two bounds
% htmp=h(find( (h>lbound) & (h<ubound)));

% determine mean & std, and amplitude of gaussian

h_ave=mean(h);
h_std=std(h);
amp=max(N);

% plot approximate gaussian fit:
hold on
xtmp=hstart:.2:hend;
plot(xtmp,amp*exp(-(xtmp-h_ave).^2/2/h_std^2),'r');
xlabel('Nanoparticle Height (nm)');
ylabel('Nanoparticle Number');
hold off

fprintf(1,'Mean (after outlier removal): %4.1f\n' ,h_ave);
fprintf(1,'Std Dev (after outlier removal): %4.1f\n' ,h_std);
