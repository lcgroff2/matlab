%SSDProb.m

%Teams A-E and Triangular Probability Inputs:

A = [7 11 16] ;
B = [1 2 6] ;
C = [4 8 12] ;
D = [4 6 11] ;
E = [4 6 8] ;
goal = 21;
teams = ['A' 'B' 'C' 'D' 'E'];
avg = zeros(size(teams));
for i = 1:length(teams);
    avg(i) = mean(teams(i));
    tot = sum(avg);
end
complete = tot<goal;
print(complete)