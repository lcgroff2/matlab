function II = loadspe(fnam)

% LOADSPE -- Load an SPE file into MATLAB

% Based on ROPER_ASCII_TO_MAT_BIN written by Igal Brener, 2003,
% Revised: Tzu Liang Loh, MIT, 2004.
% LOADSPE by Jason McNeill, Clemson, 2008.

% example:
% II = loadspe('data.spe');
% imagesc(II(:,:,1))

fid = fopen(fnam,'r'); 

if fid > 0
  header = fread(fid,2050,'uint16=>uint16'); % 2050 uint16 = 4100 bytes = 32800 bits
	
  Xdim = header(22);
  Ydim = header(329);
  Zdim = header(724);
  DataType = header(55);
  
  switch DataType
   case 0	% FLOATING POINT (4 bytes / 32 bits)
    ImMat = fread(fid,inf,'float32=>float32');
   case 1	% LONG INTEGER (4 bytes / 32 bits)
    ImMat = fread(fid,inf,'int32=>int32');
   case 2	% INTEGER (2 bytes / 16 bits)
    ImMat = fread(fid,inf,'int16=>int16');
   case 3	% UNSIGNED INTEGER (2 bytes / 16 bits)
    ImMat = fread(fid,inf,'uint16=>uint16');
  end
  
  fclose(fid);
  
  II.image = reshape(ImMat,Xdim,Ydim,Zdim);
  clear ImMat; %clear some memory
  
  %permute the X and Y dimensions so that an image looks like in Winview
  %II = permute(II,[2,1,3])
  % I commented because permute takes some time, slows it down a bit
  
else
  error('File not found');
  II=0;
end


