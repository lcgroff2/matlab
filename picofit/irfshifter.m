function yout=irfshifter(T0)
  % new irfshifter, does interpolation
  global IRFx IRFy
  dx=IRFx(2)-IRFx(1);
  shiftidx=round(T0/dx);
  lnirf=length(IRFy);
  %if shiftidx>0,
  %   yout=[zeros([1 shiftidx]) IRFy(1:(lnirf-shiftidx))];
  %elseif shiftidx<0,
  %   shiftidx=-shiftidx;
  %   yout=[IRFy((1+shiftidx):lnirf) zeros([1 shiftidx])];
  %else
  %   yout=IRFy;
  %end 
  yout=interp1(IRFx,IRFy,IRFx-T0,'*cubic',0); % extrapolation value=0
