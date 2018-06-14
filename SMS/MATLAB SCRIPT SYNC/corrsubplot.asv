% Plot data from autocorrelation and 2D, 3D trajectories
% Can run solely if the workspace data is still in use in order to plot
% more tracks.

tracknum = input('Plot tracking #');
tracknuminlist = find(goodtracks==tracknum);
mytrack = corrtracklist{1,tracknuminlist}; %Select what particle will be used

%%
% plot range for autocorrelation, 1:200 sounds good, can edit if needed
ran=1:1000;

% ranges corresponding to "early" versus "late" part of
% trajectory.  **PLEASE EDIT**  For longer trajectories, need to
% set "late=5001:6000" or similar.
%early=1:2000;
%late=2001:4000;
xx = mytrack(:,1) * pxlsize;
yy = mytrack(:,2) * pxlsize;
ww = mytrack(:,3);
ii = mytrack(:,4);
normI = ii/mean(ii);
tlen=length(xx);
tt=1:tlen;
tt=tt*rnddwelltime_s;
xx = xx-mean(xx);
yy = yy-mean(yy);
ww = ww-mean(ww);

%RMSD vs Tau
xmsd=msd(xx);
ymsd=msd(yy);
rmsd=xmsd+ymsd;

%RMSD vs Time
dx=xx(2:end)-xx(1:(end-1));
dy=yy(2:end)-yy(1:(end-1));
r_rms=sqrt(mean(dx.^2+dy.^2));
% fprintf(1,'RMSD = %7.4f pixel/frame\n',r_rms);

% plot the trajectories
figure(1);
% make a tall figure
%scrsz = get(0,'ScreenSize');
%set(gcf,'position',[10 102 900 700]);
subplot(2,2,[1 2]) % stretch first plot all the way across figure.
% first subplot is trajectory
[ax,Ilsvst] = plotyy(tt,[xx yy],tt,normI);
set(ax(1),'ylim',[-20 20],'ytick',[-20 -10 0 10 20],'box','off'); 
set(ax(2),'ylim',[0.5 1.5],'ytick',[0.5 1 1.5]); 
xlim([0 max(tt)])
set(get(ax(1),'ylabel'),'string','Position - Mean (nm)');
set(get(ax(2),'ylabel'),'string','Intensity/Mean');
xlabel('Frame number')
title(['Intensity and position trajectories of track #' num2str(tracknum)])
legend('X','Y','I')

% (2) plot autocorrelations
subplot(2,2,3)
% iacearly=acft(ii(early));
% iacearly=iacearly(2:end); % skip first point, which is just the sum
% iacearly=iacearly-min(iacearly);
% iacearly=iacearly(ran);
% iacearly=iacearly/iacearly(1); % normalize
% iaclate=acft(ii(late));
% iaclate(2:end);
% iaclate=iaclate-min(iaclate);
% iaclate=iaclate(ran);
% iaclate=iaclate/iaclate(1);
Iac=acft(ii);
Iac=Iac(2:end);
% Iac=Iac-min(Iac);
Iac=Iac(ran);
Iac=Iac-min(Iac); %Subtracting mean puts Iac on same scale as Xac, Yac.
Iac=Iac/Iac(1);
Xac=acft(xx);
Xac=Xac(2:end);
Xac=Xac-min(Xac);
Xac=Xac(ran);
Xac=Xac/Xac(1);
Yac=acft(yy);
Yac=Yac(2:end);
Yac=Yac-min(Yac);
Yac=Yac(ran);
Yac=Yac/Yac(1);
plot(ran,Xac,ran,Yac,ran,Iac)
legend('X','Y','I');
grid on
title('Autocorrelation, normalized')
ylabel('a.u.')
xlabel('Tau (frame)')

% Plot of MSD(tau)
subplot(2,2,4)
plot(ran,xmsd(ran),ran,ymsd(ran)) %,ran,rmsd(ran)/2
grid on
title('MSD(tau)')
ylabel('MSD nm^2')
xlabel('Tau (frame)')
legend('X','Y') %,'R/2'

figure(2)
%set(gcf,'position',[0 0 600 600]);
subplot(1,2,1)
twoD = plot(xx,yy,':r');
set(twoD,'Marker','.','MarkerEdgeColor','k')
axis square
grid on
title(sprintf('2D trajectory of track #%i with av RMSD %2.1f nm/frame\n', tracknum, r_rms))
xlabel('X (nm)')
ylabel('Y (nm)')
xlim([-15 15])
ylim([-15 15])

subplot(1,2,2)
threeD = plot3(xx,yy,tt,':r');
set(threeD,'Marker','.','MarkerEdgeColor','k')
box on
grid on
title(sprintf('3D trajectory of track #%i', tracknum))
xlim([-15 15])
ylim([-15 15])
xlabel('X (nm)')
ylabel('Y (nm)')
zlabel('Time (s)')
zlim([0 max(tt)])

% save data to a MATfile
fnam2 = strtrim(fnam); % remove any spaces
%fnam2=lower(fnam2); % switch to lower case
% find location of '.spe'
fnam2 = strrep(fnam2, '.SPE', ''); fnam2 = strrep(fnam2, '.sif', ''); fnam2 = strrep(fnam2, '.', '_');

saveas(figure(1),[[pathnam fnam2 '-track-' num2str(tracknum)] '-XYI-Autocorr']);
saveas(figure(2),[[pathnam fnam2 '-track-' num2str(tracknum)] '-2D+3D']);
