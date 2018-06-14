function c=myconv(a,b)

% FFT-based convolution.  No error checking, built for speed

lna=length(a);

% pad the function by at least 10%, make
% the length a power of 2.

ftlen=2^ceil(log(lna*1.3)/log(2));

aa=zeros([1 ftlen]);

bb=aa;


aa(1:lna)=a(:);
bb(1:lna)=b(:);

% hack to deal with long exponents and wraparound:
wrlen=floor((ftlen-lna)/2);
aa(lna:(lna+wrlen))=aa(lna);

%aft=fft(aa);
%bft=fft(bb);

%aft=aft(:);
%bft=bft(:);
% calculating the ifft using the inlined version:
% y = conj(fft(conj(x)))/length(x).
c=real(conj(fft(conj(fft(aa).*fft(bb)))))/ftlen;

