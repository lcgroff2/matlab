%CSV reader, TRIR Surface Plot from .csv file and SVD analysis. LCG 6/13/2016

%Cleanup from previous runs:
clear A M t cm1 i a x y U S V

%Input .csv data filename and baseline filename here:
datafile = 'lcg031017-Ge400-6750-30fs-2kshots.csv';
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

%Index out time column and wavenumber row, keeping only the data matrix.
%Changed from for loop line-by-line index method.
A = M(2:end,end:-1:2);

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
[U,S,V] = svd(A); %singular value decomposition of data matrix A into spectral components (U), Temporal components (V), and their weighted intensity (S)
V = V'; %Transpose V matrix (V --> V^T).
Snorm = diag(S)./max(diag(S)); %Normalize the Intensity components to get relative percentage of signals.

figure
plot([1:1:length(diag(S))],diag(S)./sum(diag(S)),'-o',[1:1:length(Snorm)],Snorm,'-o');
smin = min(Snorm);
smax = max(Snorm);
axis([1 5 smin-.05*smin smax+.05*smax])
legend('Normalized by Integrated Area','Normalized by Raw Amplitude')
title('First 5 Diagonal values of S Matrix')

%plot SVD components' spectra corresponding to >5% of the signal according
%to S matrix:
%To reduce the number of times we run the for loop (rather than loop over
%the whole matrix and check if a component is > 0.05), use find function to
%find the indices that are > 0.05 (usually only 2 or so).
Sidx = find(Snorm > 0.05); 
for j = 1:1:length(Sidx);
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
