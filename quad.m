function [x1, x2, flag] = quad(a, b, c)
x1 = ((-b) + sqrt((b*b-4*a*c)))/(2*a);
x2 = ((-b) - sqrt((b*b-4*a*c)))/(2*a);
if isnan(x1) && isnan(x2); %any solution and no solution yield NaNs for x1,x2
    if c ~= 0; %no solution case
        flag = 0;
    else       %any solution case
        flag = 99; 
    end
elseif isreal(x1) || isreal(x2); %is x1 or x2 real?
    if isreal(x1) && isreal(x2); % if yes, are they both real?
        flag = 2;
    else
        flag = 1;
    end
else %if neither are real, two complex roots.
    flag = 2;
end
end
