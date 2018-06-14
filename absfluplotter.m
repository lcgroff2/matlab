%Corrected Absorbance/Fluorescence data plotting script (LCG) 08/2011

color = ['b' 'g' 'r' 'c' 'm' 'y' 'k' ':b' ':g' ':r' ':c' ':m' ':y' ':k'];
i = input('how many samples do you have? ');
l = input('add legend to plot? (y/n) ','s');
filetype = input('is this absorbance (abs) or fluorescence (flu) data?:' ,'s');
c = input('Do you want corrected spectra? (y/n) ' ,'s');
figure

for N = 1:i;
    if filetype == 'abs'
        [fnam, pathnam, filterindex] = uigetfile('*.txt', 'pick your .txt file');
        cd(pathnam);
        fid = fopen(fnam,'r');
        sample = textscan(fid,'%f %f','delimiter', ',', 'headerlines',2);
        sample_x = sample{1,1};
        sample_y = sample{1,2};
        spectrum_x = sample_x;
        if c == 'y';
            spectrum_y = sample_y-min(sample_y);
        else
            spectrum_y = sample_y;
        end
        hold on
        if N <= 8
            plot(spectrum_x,spectrum_y,color(N));
        else
            plot(spectrum_x,spectrum_y,color(N:N+1));
        end
        xlabel('Wavelength (nm)');
        if c == 'y';
            ylabel('Corrected Absorbance (au)');
        else
            ylabel('Absorbance (au)');
        end
        if l == 'y';
            Labels{1,N} = input('input label for legend: ','s');
        end
    else
        if filetype == 'flu'
            [fnam, pathnam, filterindex] = uigetfile('*.txt', 'pick your .txt file');
            cd(pathnam);
            fid = fopen(fnam,'r');
            sample = textscan(fid,'%f %f','headerlines',4);
            sample_x = sample{1,1};
            sample_y = sample{1,2};
            if c == 'y';
                spectrum_y = sample_y-min(sample_y);
            else
                spectrum_y = sample_y;
            end
            hold on
            if N < 8
                plot(sample_x,sample_y,color(N));
            else
                plot(sample_x,sample_y,color(N:N+1));
            end
            xlabel('Wavelength (nm)');
            if c == 'y';
                ylabel('Corrected Fluorescence Intensity (counts/s)');
            else
                ylabel('Fluorescence Intensity (counts/s)')
            end
            if l == 'y';
                Labels{1,N} = input('input label for legend: ','s');
            end
        else
            break
        end
    end
    if N == i
        Title = input('Input title: ','s');
        title(Title);
        if l == 'y';
            legend(Labels);
        end
        hold off
    end
end