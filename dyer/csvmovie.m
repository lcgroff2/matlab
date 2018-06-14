%CSV reader, TRIR Surface Plot from .csv file and SVD analysis with autocorrelation. LCG 6/13/2016

%Cleanup from previous runs:
clear A M t cm1 i a x y U S V autoU autoV

%Input .csv file name here
datafile = 'lcg040717-Ge-t0-400nm2mW-6950nm0.7mW-100fs-1k.csv';
M = csvread(datafile);
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
    plot(cm1,A(:,i));
    xlabel('frequency (cm^-^1)')
    ylabel('Delta A')
    N(i) = getframe;
end

%Put together individual spectra as frames of movie.
for j = 1:1:length(t);
    caxis auto
    plot(cm1,A(:,j),cm1,zeros(length(cm1)),'-.k');
    axis([cm1(1) cm1(end) -1.2e-3 1.2e-3]);
    xlabel('frequency (cm^-^1)')
    ylabel('Delta A')
    title('Time Lapse of TRIR Spectra at 10 frames/s')
    N(j) = getframe;
end

%Play movie
numreps = 1; %Repeat count
fps = 20; %Frame rate
movie(N,numreps,fps)