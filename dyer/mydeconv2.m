function c=myconv2(a,b)

% FFT-based convolution.  No error checking, built for speed

lna=length(a);

% pad the function by at least 10%, make
% the length a power of 2.

ftlen=2^ceil(log(lna*1.1)/log(2));

aa=zeros([1 ftlen]);

bb=aa;

% properly pad on both sides, to deal with wraparound
lpad=floor(0.04*lna);
rpad=lna+lpad;
aa(lpad:(lna+lpad-1))=a(:);
bb(lpad:(lna+lpad-1))=b(:);
% left padding
aa(1:lpad)=a(1);
bb(1:lpad)=b(1);
% right padding
aa(rpad:end)=a(end);
bb(rpad:end)=b(end);

%aft=fft(aa);
%bft=fft(bb);

%aft=aft(:);
%bft=bft(:);
% calculating the ifft using the inlined version:
% y = conj(fft(conj(x)))/length(x).
c=real(conj(fft(conj(fft(aa)./fft(bb)))))/ftlen;

% strip off the padding
c=c(lpad:(lna+lpad-1));

%lna
%ftlen
%lpad
%rpad

%stop