%trircorrect.m Baseline correction script
Acorr = zeros(size(A(:,1)));
for i = 1:length(t);
    p = polyfit(cm1,A(:,i),2);
    Acorr(:,i) = A(:,i)-(p(1)*cm1.^2+p(2)*cm1+p(3));
end