i = 1:1:5;
a = zeros([1 10]);
for i = 1:5;
    j = 1;
    if a(j) < 5
        if i > 1
            i = i-1;
        end
        if j > 1
            a(j) = a(j-1)+1;
            j = j+1;
        else
            a(j) = a(j)+1;
            j = j+1;
        end
    end
end
a