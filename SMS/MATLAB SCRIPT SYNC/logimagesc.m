% logimage -- does a log-image, with a proper colorbar

% assumes 11-bit or 16-bit integer image data, or a sum of image frames

%function logimagesc(img)

% no zeros or negatives allowed
%img=img-min(min(img))+2;




cmax=max(max(img));
cmin=min(min(img));


% heuristic to make tick marks on color bar

norders=(log(cmax)-log(cmin))/log(10);

if norders<3
    base=2;
    % small numbers, so use base 2 scale
    lowtick = ceil(log(cmin)/log(2));
    hightick = floor(log(cmax)/log(2));
    logticks=lowtick:hightick;
    ticks=2.^logticks;
    lgimg=log(img)/log(2);
else % larger numbers, so use decimal scale
    base=10;
    lowtick = ceil(log(cmin)/log(10));
    hightick = floor(log(cmax)/log(10));
    logticks=lowtick:hightick;
    ticks=10.^logticks;
    lgimg=log(img)/log(10);
end

imagesc(lgimg)
h=colorbar; % ('YTick',logticks,'YTickLabel',ticks);
ticks=get(h,'ytick');
if base==2
  set(h,'yticklabel',2.^ticks);
else
  set(h,'yticklabel',10.^ticks);
end