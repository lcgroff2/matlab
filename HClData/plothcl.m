%plothcl: script to plot HCl CSV data from FTIR

filename = 'hcl3.csv'; %give filename of CSV file including .csv
data = csvread(filename);
xdata = data(:,[1]);
ydata = data(:,[2]);
plot(xdata,ydata)
xlabel('frequency (cm^-^1)')
ylabel('relative intensity (au)')
title('FTIR Spectrum of HCl')