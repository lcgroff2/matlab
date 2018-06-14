

theodatastrexp=theodata;
biexpresults=[];

% disp(['Running biexp fit number ' num2str(I) ])
% DataYtmp=theodatastrexp+(stdev*randn(size(theodata)));
% DataY(dran)= DataYtmp;
% montepico2
% biexpresults=[biexpresults; fitinfo.fun{1}.curval fitinfo.fun{2}.curval relintamps];

for I=1:20,
    disp(['Running biexp fit number ' num2str(I) ])
    DataYtmp=theodatastrexp+(stdev*randn(size(theodata)));
    DataY(dran)= DataYtmp;
    montepico2
    biexpresults=[biexpresults; fitinfo.fun{1}.curval fitinfo.fun{2}.curval relintamps];
end

