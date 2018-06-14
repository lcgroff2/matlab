
% corrplot.m

% Do autocorrelation analysis / plotting of trajectories.

% Does the analysis of trajectory that is in "mytrack".  Should use
% "vibcomp" first, to select "good" trajectories and do vibration
% correction.

% Plots:

% 1. A "plotyy" plot of x position, y position trajectories
% (subtracting the mean, and using a 2 by N matrix so both X and Y can
% be plotted), and intensity as the second y axis.  Should label axes
% and lines (legend) for clarity.
%
% Note: look at this plot carefully, you should discard any
% trajectories that have several "hops" larger than 0.5 pixel.  Also,
% if there is too much photobleaching (a judgment call, depends on
% the goal of the experiment), then perhaps you should truncate the
% trajectory.  There is an option for truncating the trajectory in
% "vibcomp". 

% For width trajectory:  Might give additional information, and also
% provide criteria for throwing out trajectories with strange width
% behavior (width changes by more than .5 pixel between frames,
% etc).  Note that "width" is not meaningful for image data
% that has been bandpass filtered (such as using "spetrk4").  Width
% data obtained by "spetrkg" is usually reliable.

% 2. Autocorrelation of "early" intensity dynamics and "late"
% intensity dynamics, plotted together, normalized.  Also
% position autocorrelation.

% This helps determine two things: (a) do the dynamics change over
% time (more or less blinking as photobleaching progresses), (b) Do
% the timescales of position fluctuations and intensity
% fluctuations coincide?  If so, that helps confirm the notion that
% the same underlying phenomena are involved.

% Also, it will be interesting to compare how the autocorrelation
% dynamics depend on excitation intensity.

% Note/future: Maybe should add exponential fits?

% 3. MSD plot.  Look for typical confined diffusion behavior.

% Note/future: Maybe do a fit to "confined diffusion" model?

%%
pxlsize = 65.85; %65.85nm/pixel

%% OPEN FILE
if ispc || ismac, %ISPC--True for the PC (Windows) version of MATLAB.ISPC, returns 1 for PC (Windows) versions of MATLAB and 0 otherwise.
  [fnam, pathnam, filterindex] = uigetfile('*.mat', 'Pick a MAT-file');
  %cd(pathnam); 
else
  ls *.mat
  fnam=input('enter SIF file name ','s');
end
disp('FILE')
fprintf(1,'%s\n',fnam);

load([pathnam fnam]);
[numframes,x,ntrack] = size(tracklist{1,1});
vibcomp;
if length(goodtracks)<2,
  return;
end
corrsubplot;
