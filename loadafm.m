function afm = loadafm(fnam)

% LOADAFM -- Load a Quesant AFM file into MATLAB

% Based (loosely) on quesant.c from Gwyddion source code.

% example:
% [Z,R] = loadspe('data.spe');
% imagesc(Z)

% Fileinfo structure from Gwyddion:
%typedef struct {
%    guint32 desc_offset;       // offset of description (unused)
%    guint32 date_offset;       // date location (unused)
%    guint32 palette_offset;    // unknown 
%    guint32 keys_offset;       // some sort of description again
%    guint32 image_offset;      // offset of image data (first int16 
%                               // is number of rows and cols)
%    guint32 img_p_offset;      /* offset of Z axis multiply, this points to a
%                                  relatively long block, the Z factor just
%                                  happens to be there */
%    guint32 hard_offset;       /* offset of X/Y axis width/height, this points
%                                  to a relatively long block, the factors
%                                  just happen to be there */
%    guint32 short_desc_offset; // offset of short desc (unused)

%    /* Read data */
%    guint32 img_res;           // size of image
%    gdouble real_size;         // physical size
%    gdouble z_scale;           // z-scale factor
%    const guint16 *image_data;
%    gchar *title;
%} FileInfo;

% some constants
MAGIC = [2 0 0 0 1 0 0 0];
MAGIC=MAGIC(:);

MAGIC_SIZE=length(MAGIC); % = 8
HEADER_SIZE=hex2dec('148'); % = 328 decimal


fid = fopen(fnam,'r'); 

afm=[];

if fid<=0
  disp('file not found')
  return
end

% read in the full file for debugging purposes
filedata = fread(fid,HEADER_SIZE*2,'uint8=>char');

fseek(fid,0,-1); % rewind the file

% read in the magic
mag=fread(fid,8);

if mag==MAGIC
  disp('magic OK--appears to be AFM file')
end

% header appears to be organizes as 4-byte "keys" followed by 4
% bytes containing the value

for j=1:((HEADER_SIZE-MAGIC_SIZE)/8),
  key=fread(fid,4,'uint8=>char');
  key=key(:)'; % make a row vector, else switch/case doesn't work
  value=fread(fid,1,'uint32');
  %fprintf(1,'key: %s   value: %i\n',key,value)
  switch key
   case 'SDES'
    afm.short_desc_p=value;
   case 'DESC'
    afm.desc_p = value;
   case 'DATE'
    afm.date_p = value;
   case 'PLET'
    afm.palette_p = value;
   case 'IMAG'
    afm.image_p = value;
   case 'HARD'
    afm.hard_p = value;
   case 'IMGP'
    if value>0,
      afm.image_param_p = value;
    end
  end
end

if exist('afm.short_desc_p'),
  afm.short_desc=filedata((afm.short_desc_p+1):(afm.short_desc_p+40))'; % short description (text)
else
  afm.short_desc='NO DESCRIPTION';
end
  
  
% load image resolution (number of lines?)
fseek(fid,afm.image_p,-1);
afm.image_res=fread(fid,1,'uint16');

% load image "real" size (microns?)
fseek(fid,afm.hard_p,-1);
afm.real_size=fread(fid,1,'float'); % 32-bit float

% load value scale factor (?)
fseek(fid,afm.image_param_p+8,-1);
afm.z_scale = fread(fid,1,'float',4);

% load 16-bit data
fseek(fid,afm.image_p+2,-1); % if I am reading Gwyddion code correctly, 
% raw data is right after image_p
afm.raw_data=fread(fid,afm.image_res^2,'uint16');
afm.zdata=reshape(afm.raw_data,afm.image_res,afm.image_res)'*afm.z_scale;
afm.zdata=fliplr(afm.zdata)*1000; % flip to match Gwyddion, convert to nm

afm.xy_range=linspace(0,afm.real_size,afm.image_res);

%


fclose(fid);
return


