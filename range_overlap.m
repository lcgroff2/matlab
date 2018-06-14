function logic_result = range_overlap(x)
if length(x) ~= 5;
    error('Invalid Length. vector must contain 5 elements (min1, max1, min 2, max2, overlap range)');
else
    %store values of vector as appropriate variables
    min1 = x(1);
    max1 = x(2);
    min2 = x(3);
    max2 = x(4);
    range = x(5);
end
width1 = max1 - min1;
width2 = max2 - min2;
max2min1 = max2 - min1;

if (min2 > max1) || (min1 > max2); % no overlap condition, 
    logic_result = 0;
elseif max2 > max1; %if overlap is between max1 and min2, this yields true
    overlap = max1-min2;
    if overlap >= (range)-1;
        logic_result = 1;
    else
        logic_result = 0;
    end
else %overlap lies between max2 and min1
    overlap = max2-min1;
    if overlap >= (range)-1;
        logic_result = 1;
    else
        logic_result = 0;
    end
end
end