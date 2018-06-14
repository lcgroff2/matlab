%signal averaging example:
x = 0:1:128;
gaus = jmgauss(x-0.5*max(x),0.5*max(x));
sig = jmgauss(x-0.5*max(x),0.5*max(x))-0.5*max(gaus)*sin((x-0.5*max(x))/15);

noigaus = zeros(size(gaus));
noisig = noigaus;
pershot = noigaus;
summedgaus = noigaus;
summedsig = noigaus;
sumpershot = noigaus;
randn('state', sum(100*clock));

for i = 1:2000;
    noise = 0.1*max(gaus)*randn(size(gaus))+0.5;
    noigaus(i,:) = gaus+noise;
    noisig(i,:) = sig+noise;
    pershot(i,:) = -log10((noisig(i,:)./noigaus(i,:)))-(-log10((noigaus(i,:)./noigaus(i,:))));
    summedgaus = summedgaus+noigaus(i,:);
    summedsig = summedsig+noisig(i,:);
    sumpershot = sumpershot+pershot(i,:);
end
summedshots = -log10(((summedsig./summedgaus)))-(-log10((summedgaus./summedgaus)));

figure
subplot(4,1,1);
plot(x,noigaus(i,:),x,noisig(i,:));
subplot(4,1,2);
plot(x,-log10((noisig(i,:)./noigaus(i,:)))-(-log10((noigaus(i,:)./noigaus(i,:)))));
subplot(4,1,3);
plot(x,sumpershot./2000);
subplot(4,1,4);
plot(x,summedshots./2000);
clear gaus sig noise pershot noigaus noisig