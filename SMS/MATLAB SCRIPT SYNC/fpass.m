
function imfilt=fpass(img,H)

ftimg=fft2(img);
ftimg=(fftshift(ftimg)'); % bring low freq to the center
% multiply by filter, shift back, and ifft
imfilt=ifft2(ifftshift(H.*ftimg));


imfilt=rot90((real(imfilt))',2); %imfilt=fliplr(flipud((real(imfilt))'));

% kludge for off-by-one error

%imfilt(1:(end-1),1:(end-1))=imfilt(2:end,2:end);
imfilt(2:end,2:end)=imfilt(1:(end-1),1:(end-1));
imfilt=imfilt-min(min(imfilt));