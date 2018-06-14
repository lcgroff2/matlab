function y=loadpico(fnam)
fid=fopen(fnam,'r');

y=fscanf(fid,'%f',[1,inf]);

fclose(fid);

