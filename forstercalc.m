%Forster Radius Calculator Using PhotoChemCAD Spectral Data
%04/16/12 - LCG

%Define constants for the R0 calculation:

Qd = 0.66; %Donor quantum yield (unitless)
n = 1.90; %Refractive index of medium (unitless)
k2 = 2/3; %Transition dipole orientation factor (assume 2/3)
dl = 1; % spacing between data points in spectrum.
conc = 2.0965e-8; %Molar concentration of acceptor (if Shimadzu data)
ovrl = 460; % Spectral overlap lower bound (nm)
ovru = 540; % Spectral overlap upper bound (nm)
vis = 'shim'; %type 'shim' for Shimadzu UV-Vis data, 'cad' for PhotoChemCAD
              %data
flu = 'pti';  %type 'pti' for PTI fluorescence data, 'cad' for PhotoChemCAD
              %data

wav = ovrl:dl:ovru; %Wavelength range of spectral overlap
wav = wav';

%Give spectra files:
%Donor emission spectrum:
don = 'lcg051014-0%THF-0.058-450-0.50-473.txt';
%Acceptor UV-Vis Spectrum (abs units) or PhotoChemCAD extinction spectrum
%(M-1*cm-1 units):
acc = 'lcg051014-0%THF-PFBTCPNs-450.txt';

fid = fopen(don,'r');
CAD = strcmp(flu,'cad');
PTI = strcmp(flu,'pti');
if PTI == 1
    donor = textscan(fid,'%f %f','headerlines',4);
else
    if CAD == 1
        donor = textscan(fid,'%f %f','headerlines',16);
    else
        disp('invalid fluorescence data type');
    end
end
donorx = donor{1,1};
donory = donor{1,2};
donory = donory-min(donory); %Corrected donor emission spectrum needed.
fclose(fid);

%Truncate donor emission spectrum to overlap range:
donory = donory(find(donorx==ovrl):find(donorx==ovru));

fid = fopen(acc,'r');
SHIM = strcmp(vis,'shim');
VCAD = strcmp(vis,'cad');
if SHIM == 1
    acceptor = textscan(fid,'%f %f','delimiter',',','headerlines',2);
else
    if VCAD == 1
        acceptor = textscan(fid,'%f %f','headerlines',23);
    else
        disp('invalid absorbance data type');
    end
end
acceptorx = acceptor{1,1};
acceptory = acceptor{1,2};
if VCAD == 1
    acceptorx = acceptorx(2:2:end); %data spacing is originally 0.25 nm,
    acceptory = acceptory(2:2:end); % adjust to dl = 0.5 nm
    acceptorx(1:2:end) = round(acceptorx(1:2:end));%round whole wavelengths
    % from 0.99 to integers.
else
    if SHIM == 1
%         if (acceptorx(2)-acceptorx(1)) < dl %In case PTI data is 1 nm step
%             acceptorx = acceptorx(1:2:end); % and Shimadzu data is 0.5 nm
%             acceptory = acceptory(1:2:end); % step, this corrects it
%         end
        acceptory = acceptory/conc; %Convert to extinction spectrum,
        %assume 1 cm path length (M-1*cm-1)
    end
end
fclose(fid);

% Truncate acceptor extinction spectrum to overlap range:
acceptory = acceptory(find(acceptorx==ovrl):find(acceptorx==ovru));

%Plot spectra to check:
figure;
plot(wav,donory,wav,acceptory);
disp('Does the overlap look correct? If so, press return else hit ctrl+c');
pause

%J integral calculation (units of M-1 cm-1 nm^4):
J = sum(donory.*acceptory.*wav.^4*dl)/sum(donory*dl)

%Lakowicz Equation 13.5 (p. 466) yields R0 in Angstroms:
R0 = 0.211*(k2*n^-4*Qd*J)^(1/6)