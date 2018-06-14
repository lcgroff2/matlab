
% sif2mat - loads each sif file in a folder, and writes out as a matfile of the
% same name, but with the .mat extension. The resulting .mat file is smaller by a factor
% of 2 (thus easier to copy from server, etc), and loads ~5x faster than the sif file

% It only does the conversion if the matfile does not already exist.

% note 1: this script does not erase or delete the original sif file
% note 2: matfile is smaller, because it saves image data as 16 bit unsigned integer (uint16),
% whereas the Solis software saves the data as 32-bit floating point (dumb!)

names=dir('*.sif');

for idx=1:length(names)
    sifnam=names(idx).name;
    matnam=sifnam;
    matnam(end-2:end)='mat';
    if exist(matnam,'file')
        fprintf(1,'%s exists, skipping...\n',matnam);
    else
        fprintf(1,'%s => %s\n',sifnam,matnam);
        data=sifread(sifnam);
        save(matnam,'data');
    end
end
