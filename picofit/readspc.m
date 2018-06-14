function outdata = ReadSpc(file)
if ~isstr(file)
   error('readspc: Filename must be a string');
end
f = fopen(file,'r');
if f == (-1)
   error('readspc: Cound not open file');
end
fread (f,768,'char');
outdata = fread(f,2048,'float32');
%   outdata = flipud(outdata);
fclose(f);
