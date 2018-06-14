function da1 = despike(da)
% this function get rid of spikes in transient
% da should be one dimentional array

dev=sqrt(mean(da))*5;
if dev < 30
   dev= 30;
end
da1(1)=da(1);
da1(length(da))=da(length(da));
for i = 2:length(da)-1
   md=(da(i-1)+da(i+1))/2;
   if (da(i) - md) > dev
      da1(i) = md;
   else
      da1(i) = da(i);
   end
end
