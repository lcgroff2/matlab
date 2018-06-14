mwH2O = 18; %g/mol
mwTHF = 72; % g/mol
PvapTHF = 143; %vapor pressure of pure THF at 298 K (torr)
PvapH2O = 20; %vapor pressure of pure H2O at 298 K (torr)
R = 62.3; %L torr mol-1 K-1
T = 298; %K
V = 1; %L

%initial mass of THF/H2O assuming 10 mL 10% w/w THF in H2O ~10 g solution:
mTHF = 1;
mH2O = 9;
concTHF = (mTHF/(mTHF+mH2O))*100;
N = 0;
test = concTHF <= 1;

while test == 0;
    N = N+1;
    %initial moles (ni) and mole fractions in liquid (x):
    niTHF = mTHF/mwTHF;
    niH2O = mH2O/mwH2O;
    nitotal = niTHF+niH2O;
    xTHF = niTHF/nitotal;
    xH2O = niH2O/nitotal;
    
    %Partial Pressures in vapor (Raoult's Law) and mole fractions in vapor (y):
    pTHF = PvapTHF*xTHF;
    pH2O = PvapH2O*xH2O;
    ptotal = pTHF+pH2O;
    yTHF = pTHF/ptotal;
    yH2O = pH2O/ptotal;
    
    %Assume ideal gas (PV = nRT):
    nTHF = (pTHF*V)/(R*T);
    nH2O = (pH2O*V)/(R*T);
    mTHF = mTHF-nTHF*mwTHF;
    mH2O = mH2O-nH2O*mwH2O;
    
    concTHF = (mTHF/(mTHF+mH2O))*100;
    test = concTHF <= 1;
end
Extractions = N