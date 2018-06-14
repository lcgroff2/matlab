
% fadplot.m -- plots fluorescence anisotropy decay data.

% JDM 5/2013

% NOTE: THIS CODE DOES NOT FIT DATA. It just plots the r(t) based on the data.


% User must edit the following lines.

% zero degree (parallel) polarization
DataFile0='lcgjyf073013-40%-para-s2.asc';
%IRFFile0='lcg062113-ud-para-i1.asc'; 
A55_0 = 334560; %KWW amplitude contribution from 0° Data

% 55 degree (magic angle) polarization
DataFile55='lcgjyf073013-40%-mag-s2.asc';
IRFFile55='lcgjyf073013-40%-mag-i2.asc'; % for now, we just use 55 degree IRF
A55 = 368084; % Amplitude from KWW fitting to magic angle trace = A55_KWW

% 90 degree (perpendicular) polarization
DataFile90='lcgjyf073013-40%-perp-s2.asc';
%IRFFile90='lcg062113-ud-perp-i1.asc';
A55_90 = 362272; %KWW amplitude contribution from 90° Data

dt=2.6333;                % dwell time in ps (picoseconds per
                          % bin). This should be re-measured
                          % periodically.

dran=1500:2500; % Data range.  Carefully edit to include 100-300 ps
                % before t0, and to not go much past 4*tau. 
ReverseMode=1; % set to "1" for reverse-mode data (PIN sent to STOP input on TAC)

% load datafiles
DataY0=loadpico(DataFile0);
DataY55=loadpico(DataFile55);
DataY90=loadpico(DataFile90);
%DataY=loadpico(DataFile);
DataX=(1:length(DataY0))*dt;

IRFy=loadpico(IRFFile55); % for now, just using 55 degree
IRFx=DataX;

if ReverseMode==1
    DataY0=DataY0(end:-1:1); % flip it
    DataY55=DataY55(end:-1:1); % flip it
    DataY90=DataY90(end:-1:1); % flip it
    IRFy=IRFy(end:-1:1); % flip it
end

% OK, now let's calculate and plot r(t)
DataX=DataX(dran);
DataY0=DataY0(dran);
DataY55=DataY55(dran);
DataY90=DataY90(dran);
IRFx=IRFx(dran);
IRFy=IRFy(dran);

%Amplitude correction factors:
A90 = A55/A55_90;
A0 = A55/A55_0;
DataY0 = A0*DataY0;
DataY90 = A90*DataY90;

% find t0
[IRFmax,IRFidx]=max(IRFy);
DataX=DataX-DataX(IRFidx);

% subtract "dark counts" and left-over intensity from previous pulse
DataY0=DataY0-min(DataY0(1:200))+1; % estimate as min of first 200 points
% DataY55=DataY55-min(DataY55(1:200))+5; % slight offset
DataY90=DataY90-min(DataY90(1:200))+1;

% find the leading edges, to try to deal with t0 drift or systematic differences in t0 at different angles
[val55,maxidx55]=max(DataY55);
[val,edgeidx55]=min((DataY55(1:maxidx55)-val55/2).^2); % leading edge as half-max
[val90,maxidx90]=max(DataY90);
[val,edgeidx90]=min((DataY90(1:maxidx90)-val90/2).^2); % leading edge as half-max
[val0,maxidx0]=max(DataY0);
[val,edgeidx0]=min((DataY0(1:maxidx0)-val0/2).^2); % leading edge as half-max

% % do the data shifting to line up the leading edges to the one at magic angle
diff0=edgeidx0-edgeidx55;
if diff0>0 % if the parallel data is shifted to the right, then shift by diff0 to the left
    tmp=DataY0((diff0+1):end); % it might be diff0+1 instead?
    DataY0(1:length(tmp))=tmp;
elseif diff0<0 % then move by -diff0 to the right
    tmp=DataY0(1:(end+diff0));
    DataY0((-diff0+1):end)=tmp;
else
    disp('no shift for parallel curve');
end
diff90=edgeidx90-edgeidx55;
if diff90>0 % if the parallel data is shifted to the right, then shift by diff0 to the left
    tmp=DataY90((diff90+1):end); % it might be diff0+1 instead?
    DataY90(1:length(tmp))=tmp;
elseif diff0<0 % then move by -diff0 to the right
    tmp=DataY90(1:(end+diff90));
    DataY90((-diff90+1):end)=tmp;
else
    disp('no shift for perp curve');
end

% %arbitrarily set DataY55 for all times earlier than slightly before the leading edge, to remove odd fluctuations
% %at early times
% DataY55(1:(edgeidx55-30))=max(DataY55)/100;

%tailmatching 0 and 90 degree data
tmatch = 500; %position in time where we want the datasets to overlap (picoseconds)
tidx = min(find(DataX > tmatch));
deltaI = DataY0(tidx)-DataY90(tidx); %Difference in intensity at tailmatching time index
DataY0 = DataY0-deltaI; %Shift parallel data by intensity difference

%sanity checking to make sure that the data matches at tidx
% figure
% plot(DataX(tidx-10:1:tidx+10),DataY0(tidx-10:1:tidx+10),DataX(tidx-10:1:tidx+10),DataY90(tidx-10:1:tidx+10))
DataY0 = DataY0(:);
DataY90 = DataY90(:);
A_diff = DataY0\horzcat(DataY90,ones(size(DataY90))); %concatenate a column of ones to represent offset
D_match = DataY0*A_diff(1)-DataY90;

% arbitrarily use last 200 points to estimate the G or alpha parameter
% g=mean(DataY0((end-200):end))/mean(DataY90((end-200):end));
g = A_diff(2);
r = D_match./(DataY0+2*g*DataY90);
r2 = D_match./(3*DataY55');
                                                      
% also added 0.005*max to S_t or DataY55, since this reduces noise before t0 and has negligible effect on the rest of curve

figure(1)
plot(DataX,r,DataX,r2)
axis tight
xlabel('Time (ps)')
ylabel('Anisotropy')
title('mldivide approach')
legend('0+2*g*90 denom.','3x magic angle denom.')

figure(2)
plot(DataX,IRFy/max(IRFy)*max(DataY55),DataX,A0*DataY0,DataX,DataY55,DataX,A90*DataY90)
axis tight
legend('IRF','0 degrees','55 degrees','90 degrees')
xlabel('Time (ps)')
ylabel('Corrected Intensity')


return

