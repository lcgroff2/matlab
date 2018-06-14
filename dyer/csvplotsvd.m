%CSV reader, TRIR Surface Plot from .csv file and SVD analysis. LCG 6/13/2016

%Cleanup from previous runs:
clear A M t cm1 i a x y U S V

%Input .csv data filename and baseline filename here:
datafile = 'lcg040518-HeminTest2k-100fs-2.csv';
M = csvread(datafile);
% basefile = 'lcg011218-hemin-test-7100nm400nm-2k-bg.csv';
% Mb = csvread(basefile);
% Mb = Mb'; %switch rows/columns so that each row is a transient and each column is a spectrum at a given time.

M = M'; %switch rows/columns so that each row is a transient and each column is a spectrum at a given time.
t = M(1,:); %first column is time in ps
t = t(2:end); %First number is an arbitrary zero, so we eliminate it.
t = t(end:-1:1); %reverse time axis so it goes low to high from left to right

cm1 = M(:,1); %IR frequencies
cm1 = cm1(2:end); %First number is an arbitrary zero, so we eliminate it.

%Preallocate empty matrix of zeros for A = U*S*V' SVD analysis.
A = zeros(length(cm1),length(t));

for i = 1:1:length(cm1);
    a = M(i+1,:); %copy current row as a vector.
    a = a(2:end); %first element of each row vector is a frequency in cm-1, eliminate it.
    a = a(end:-1:1); %flip row to match time axis
    A(i,:) = a; %store in A matrix
    
%     %Subtract Baseline 
%     b = Mb(i+1,:);
%     b = b(2:end);
%     A(i,:) = A(i,:)-mean(b); %Take average of baseline at frequency i and
%                              %subtract it from the signal at frequency i
end

%Show data matrix as 3D surface plot:
figure
imagesc(cm1,t,A') %Surface plot of data
% surf(t',cm1',A) %for short time vectors, comment/uncomment as needed
axis tight
zlabel('Delta A')
xlabel('Frequency (cm^-^1)')
ylabel('Absolute Delay Time (ps)')
title('Raw Data')
h = colorbar;
set(get(h,'title'),'string','Delta Abs (OD)','Rotation',270.0);

%Singular Value Decomposition analysis:
[U,S,V] = svd(A); %singular value decomposition of data matrix A
V = V';
Snorm = diag(S)./max(diag(S));

figure
plot([1:1:length(diag(S))],diag(S)./sum(diag(S)),'-o',[1:1:length(Snorm)],Snorm,'-o');
smin = min(Snorm);
smax = max(Snorm);
axis([1 5 smin-.05*smin smax+.05*smax])
legend('Normalized by Integrated Area','Normalized by Raw Amplitude')
title('First 5 Diagonal values of S Matrix')

%plot SVD components' spectra corresponding to >10% of the signal
for j = 1:1:4;
    figure
    subplot(2,1,1);
    plot(cm1,U(:,j),cm1,zeros(length(U(:,j))),'-.k')
    title('Spectrum of Component from Matrix U')
    xlabel('Frequency cm^-^1')
    ylabel('Delta A')
    axis tight
    subplot(2,1,2);
    plot(t,V(j,:),t,zeros(length(V(j,:))),'-.k')
%     plot(t,V(j,:),t,zeros(length(V(j,:))),'-.k') %for short time vectors,
%     comment/uncomment as needed
    title('Transient of Component from V^T')
    xlabel('time (ps)')
    ylabel('Delta A')
    axis tight
end
