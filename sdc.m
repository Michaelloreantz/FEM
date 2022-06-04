%% EERI SDC 2022
% Michael Loreantz Steven Tambunan
%% Prerequisites
clear
close all
clc
%% Directories
GM1_acc='/Users/Michael/Documents/GM1_2022_acc.csv';
GM1_spectrum='/Users/Michael/Documents/GM1_2022_Spectra.csv';
GM2_acc='/Users/Michael/Documents/GM2_2022_acc.csv';
GM2_spectrum='/Users/Michael/Documents/GM2_2022_Spectra.csv';
%% Ground Motion 1
% Time History
TH1=importdata(GM1_acc);
xTH1=TH1.data(:,1);                   %Time Values
yTH1=TH1.data(:,2);                   %Acceleration Values
fs=1/(TH1.data(3,1)-TH1.data(2,1));   %Sampling frequency, number of data in 1 second

figure(1);
subplot(3,1,1);
plot(xTH1,yTH1), title('Ground Motion Data (Time Domain)'), xlabel('Time (s)'), ylabel('Acceleration (g)'), xticks(0:2:length(xTH1)), grid;                        %Ground motion plotting
hold on;

LyTH1=length(yTH1);                     %Data length
LyTH1p2=2^nextpow2(LyTH1);               %Data length in power of 2
ff1=fft(yTH1,LyTH1p2);                    %Fast Fourier Transform
fff1=ff1(1:LyTH1p2/2);                     %Single Sided
fff1=fff1/max(fff1);                       %Normalized Amplitude
xfft1=fs*(0:LyTH1p2/2-1)/LyTH1p2;          %Frequency Axis
subplot(3,1,2);
semilogx(xfft1,abs(fff1)),title('Normalized Fourier Transform (Frequency Domain)'), xlabel('Frequency (Hz)'), ylabel('Normalized Amplitude'), grid;   %Normalized FFT Plotting
hold on;

% Spectrum Response
Spectrum1=importdata(GM1_spectrum);
xSpectrum1=Spectrum1.data(:,1);         %Period values
ySpectrum1=Spectrum1.data(:,2);         %Sa values (g)
subplot(3,1,3);
plot(xSpectrum1,ySpectrum1), title('Spectrum Response'), xlabel('Period (s)'), ylabel('Sa(g)'), xticks(0:1:length(xSpectrum1)),yticks(0:0.1:length(ySpectrum1)),grid;            %Spectrum Response Plotting
hold on;
%% Ground Motion 2
% Time History
TH2=importdata(GM2_acc);
xTH2=TH2.data(:,1);                   %Time Values
yTH2=TH2.data(:,2);                   %Acceleration Values
fs=1/(TH2.data(3,1)-TH2.data(2,1));   %Sampling frequency, number of data in 1 second

figure(2);
subplot(3,1,1);
plot(xTH2,yTH2), title('Ground Motion Data (Time Domain)'), xlabel('Time (s)'), ylabel('Acceleration (g)'), xticks(0:2:length(xTH2)), grid;                        %Ground motion plotting
hold on;

LyTH2=length(yTH2);                     %Data length
LyTH2p2=2^nextpow2(LyTH2);              %Data length in power of 2
ff2=fft(yTH2,LyTH2p2);                   %Fast Fourier Transform
fff2=ff2(1:LyTH2p2/2);                    %Single Sided
fff2=fff2/max(fff2);                       %Normalized Amplitude
xfft2=fs*(0:LyTH2p2/2-1)/LyTH2p2;        %Frequency Axis
subplot(3,1,2);
semilogx(xfft2,abs(fff2)),title('Normalized Fourier Transform (Frequency Domain)'), xlabel('Frequency (Hz)'), ylabel('Normalized Amplitude'), grid;   %Normalized FFT Plotting
hold on;

% Spectrum Response
Spectrum2=importdata(GM2_spectrum);
xSpectrum2=Spectrum2.data(:,1);         %Period values
ySpectrum2=Spectrum2.data(:,2);         %Sa values (g)
subplot(3,1,3);
plot(xSpectrum2,ySpectrum2), title('Spectrum Response'), xlabel('Period (s)'), ylabel('Sa(g)'), xticks(0:1:length(xSpectrum2)),yticks(0:1:length(ySpectrum2)),grid;        %Spectrum Response Plotting
hold on;