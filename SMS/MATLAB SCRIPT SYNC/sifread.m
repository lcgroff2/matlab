
function im = sifread(fnam)

% sifread function

% Based on information gleaned from Lautenegger's sifreader script
% but rewritten from scratch to clean it up and to handle Neo sCMOS
% data.

% usage:
% >> im=sifread('file.sif');
% >> imagesc(im.image(:,:,1))
% >> axis image

% note: you may need to rotate the plot (fliplr, flipud, rot90) or 
% the transpose operator ', or some combination, to make it match 
% up with Solis image.  For example:
% >> imagesc(im.image(:,:,1)') % flips x,y using the transpose operator '
% >> imagesc( flipud(im.image(:,:,1)) ); % flips upside down
% >> imagesc( fliplr(im.image(:,:,1)) ); % flips left-right

% note: you also should run "axis image" to get square pixels

% 2012 Jason D. McNeill

% To-Do: 
% (no known issues)

fid=fopen(fnam,'r');

if fid<0
  error('file not found')
end

hdr{1} = fgetl(fid);
hdr{2} = fgetl(fid); % header numbers: 65566 1
hdr{3} = fgetl(fid); % settings lines (includes some zeros)
hdr{4} = fgetl(fid); % detector type
hdr{5} = fgetl(fid);  % image area (and a 60 for no apparent reason)
hdr{6} = fgetl(fid); % full file name and path
hdr{7} = fgetl(fid); % 65538 2048

im.filename=hdr{6};

% 165h to 965h seem to be junk
% 966h - "65538", then more gibberish, including SR303i and SR163\r65538
% "Pixel number6\r"
% "Counts12\r"
% A72: "Pixel number65541"

% read lines up until "Pixel number"

curline='';
while isempty(strfind(curline,'Pixel number')),
  curline=fgetl(fid);
end

skip=fgetl(fid); % gives "Counts12"

hdr{9} = fgetl(fid); % image size settings
hdr{10} = fgetl(fid); % frame settings

imagesizedat = sscanf(hdr{9},'Pixel number%d %f %f %f %f %f %f %f %f');
framesizedat = sscanf(hdr{10},'%f %f %f %f %f %f %f');

im.pixperframe = imagesizedat(9);
im.totalpixels = imagesizedat(8);
im.datasize = imagesizedat([4 3 6]);
im.width = imagesizedat(4);
im.height = imagesizedat(3); % need to check rotation using a rectangular dataset
im.nframes = imagesizedat(6);

disp('REAL DATA')
fprintf(1,' Real size: %ix%ix%i\n',im.width,im.height,im.nframes);

% note: I am ignoring binning for now -- will revisit after I
% generate some test data with different x and y bin sizes.
% Maybe framesizedat(6:7) are the bin sizes
im.vbin=framesizedat(6);
im.hbin=framesizedat(7);

% next read microsecond timestamps

im.timestamps_microsec = fscanf(fid,'%d\n', im.nframes);

skipline=fgetl(fid); % some junk?
im.stamps_misc = fscanf(fid,'%d\n',im.nframes); % seem to be another set of timestamps using another unit

% first we will read into a linear array, and then do some checking
% on the first and last data points, as a check. After we have thoroughly 
% tested the script, we should combine the fread and reshape lines,
% and remove the checks.

% note that if we are off by one byte, the 32-bit floating point values will be junk, and it might be necessary to read in a byte or three to compensate, using skipbytes=fread(fid,numbytes,'uchar'); , where numbytes is the number of bytes to read. Adjust until resulting images look normal.

% might want to read into a 16-bit unsigned integer array, to save
% memory. To do this, change to 'float32=>uint16'

im.hsize = im.width/im.hbin;
im.vsize = im.height/im.vbin; % take binning into account

if im.hbin>1,
  fprintf(1,' Binning: %ix%i\n',im.hbin,im.vbin);
  fprintf(1,' Binned data size: %ix%ix%i\n',im.hsize,im.vsize,im.nframes);
else
  fprintf(1,' No binning\n');
end

imraw=fread(fid,im.hsize*im.vsize*im.nframes,'float32=>uint16');

if (imraw(1)<1) || (imraw(1)>65535)
  disp('Warning: offset problem, adjust script')
end
if (imraw(end)<1) || (imraw(end)>65535)
  disp('Warning: offset problem, adjust script')
end

% might be height, width instead - need to check
im.image = reshape(imraw,[im.hsize im.vsize im.nframes]);

%im.image = reshape(imraw,[im.height im.width im.nframes]);

%%
if length(im.timestamps_microsec)>5,
  im.dwelltime_s=(im.timestamps_microsec(5)-im.timestamps_microsec(1))/4/1e6;
  im.xdat=((1:im.nframes)-1)*im.dwelltime_s;
  fprintf(1,' Dwell time per frame: %8.6f sec\n',im.dwelltime_s);
end

framefivesec = round(5/im.dwelltime_s);
ypeakim = max(max(im.image(:,:,1:framefivesec)));
[im.maxI,im.maxIloc] = max(ypeakim);

[im.xmaxI,im.ymaxI] = find(im.image(:,:,im.maxIloc) == im.maxI);

if numel(im.xmaxI) > 1 || numel(im.ymaxI) >1
    im.xmaxI = im.xmaxI(1,1);
    im.ymaxI = im.ymaxI(1,1);
end

fclose(fid);
