
% msdsim - simulate the effect of per-frame uncertainty, number of frames,
% and diffusion coefficient (or per-frame-msd) on the RMSD curve.

% The point is to:
% (a) demonstrate that having more frames at least partially compensates for having poorer
% per-frame resolution at higher frame rates, I believe the overall tracking
% uncertanty, msd(1)-msd(0) mostly depends on the total number of photons collected in all
% frames, assuming that the other sources of noise (readout noise, etc) are minimal
% (b) familiarize students with how MSD curves vary from run to run
% (c) help test proposed experimental parameters 
% (d) provide proof to skeptics that we can examine very short "jumps" under most
% circumstances, either in a statistical way (if there are a lot of jumps occurring on fast time scales)
% (e) to help us make/refine/test wavelet analysis, etc. 

% To-do: make a "hopsim" code that simulates jumping on various length scales and timescales,
% perhaps following known transport models or ate least mimicking them to some degree. One
% possible outcome of that is to see to what extent acorr can help find length scales and time
% scales of jumps. Also to test whether wavelets or segmented autocorrelation can help analyze this
% kind of data
% Also: make a generalized hopping simulation code with various sites, positions, and interconversion
% rates, and then work on analysis code
% also test to see if "rolling" trajectories can be obtained for lots of tiny jumps, and under what
% circumstances rolling trajectories start to "look like" discrete jumps (i.e., when there are
% bottlenecks that separate two regions.
% also under what circumstances mutliple rate constants in the model can lead to one observed rate.
% also to see if "picket fence" or "island" or narrow bottlenecks can be differentiated


rmsd1 = 0.2; % random walk displacement of ~0.1 nm on average, per frame

rmsd0 = 2; % 2 nm uncertainty per frame

nframes=500;

xexact=cumsum(randn([1 nframes]))*rmsd1; % exact random walk trajectory

xnoisy=xexact+rmsd0*randn([1 nframes]); % random walk trajectory with random noise added

xmsd=msd(xexact);

noisymsd=msd(xnoisy);

ran=1:100;

plot(ran,xmsd(ran),ran,noisymsd(ran))
