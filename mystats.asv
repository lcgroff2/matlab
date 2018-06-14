function [x_bar, stdev] = mystats(x)

stdev = 0;
x_noNaNs = x; %copy x vector to leave original unaltered
nanidx = find(isnan(x) == 1); %Find indices where x(i) == NaN
x_noNaNs(nanidx) = []; %set NaN indices to empty values to remove from vector

x_bar = sum(x_noNaNs)/length(x_noNaNs); %easily calculate the mean

%piecewise calculate the standard deviation:
stdev = (x_noNaNs - x_bar).^2; %(x(i)-mean)^2
stdev = sum(stdev);            % sum the numerator
stdev = stdev/(length(x_noNaNs)-1); % divide by N-1 elements to get variance
stdev = sqrt(stdev);                %square root the variance to get stdev
end