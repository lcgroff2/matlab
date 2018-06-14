%Evanescent wave penetration depth by ATR with Ge crystal and 45 degree
%geometry for water sample:

freq = 1566:1:1730;
wav = (1./freq);
nGe = 4;
nH2O = 1.26;
theta = 45;

for i = 1:length(wav);
    d_p(i) = wav(i)./(2*pi()*sqrt((nGe*nGe*sind(theta)*sind(theta))-(nH2O*nH2O)));
end
