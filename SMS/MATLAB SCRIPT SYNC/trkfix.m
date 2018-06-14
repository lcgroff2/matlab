
% trkfixrun.m

% script for analyzing several data files sequentially

% for "fixed" particles.  This script should be more sensitive to
% dim particles.

% path to data.  Uncomment to set the directory if desired
%cd('c:\data\Tracking');

% file prefix.  Analyzes all SPE files starting with these characters
%fileprefix='jbyuFN08';

close all
clear
clc

% tracking parameters
% bandpass parameters
bp1 = 1;
bp2 = 30;
% threshold calculating parameter
threshfrac = 0.02;% a value of 0.1 means only include peaks that are 10% of the max, or higher
maxnumparticles = 8;
% minimum separation between peaks for pkfnd (peak finder) function
minsep = 6;
% number of pixels (on both sides of center) to use in gaussian fitting
gpixels=6;

showimages=0;

% end of user-adjustable parameters

cdpwd = pwd; addpath(pwd);
disp('Choose a folder.')
fdname = uigetdir('.'); disp(fdname); cd(fdname)

fstructSPE = dir('*.SPE');
fstructSIF = dir('*.sif');
fstruct = fstructSPE;
fstruct((length(fstructSPE)+1):(length(fstructSPE)+length(fstructSIF)),1) = fstructSIF;

if length(fstruct)<1,
  disp('-Bad file name.');
  return;
end


fprintf(1,'The files that will be analyzed are:\n');
for cnt=1:length(fstruct),
  fprintf(1,'%s\n',fstruct(cnt).name);
end

anfiles = input('Analyze these files? (Y/n) ','s');
if anfiles == 'n' | anfiles == 'N';
    disp('-No.')
    cd(cdpwd)
    return
elseif isempty(anfiles) | anfiles == 'y' | anfiles == 'Y';
    disp('-Yes.')
    cropsec = input('Crop frame (sec) = ');
    for cnt=1:length(fstruct),
        fnam=fstruct(cnt).name;
        trkfix2
    end
else
    disp('-Invalid answer.')
    cd(cdpwd)
    return
end

cd(cdpwd)
