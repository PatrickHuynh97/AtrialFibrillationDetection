%% Uses MIT-BIH AF Database ECG signals and tests several Baseline Wander removal methods
clear all
addpath('mcode');

% Downloading ECG Signal from Atrial Fib. Database
%[signal,Fs,tm]=rdsamp('afdb/04015',[], 1000); % # AFIB ECG
%[signal,Fs,tm]=rdsamp('nsrdb/16265',[], 10000); % # NSR ECG

[signal,Fs,tm]=rdsamp('afdb/08455',[], 1000); % # NSR ECG

ecg = signal(:,2); % select lead

figure('Name','Removal of Baseline Wander with different Algorithms','NumberTitle','off');

% Original Signal
subplot(3,1,1)
plot(tm,ecg), grid
title(['\fontsize{16}' 'Raw ECG Signal']), 
xlabel Sample, ylabel 'Voltage (mV)'

%{
%% Non-linear de-trending with Polyfit
% Code from mathworks.com/help/signal/ug/remove-trends-from-data.html
% Does not work as expected. 
opol = 6;
[p,s,mu] = polyfit(tm,ecg,opol);
f_y = polyval(p,tm,[],mu);
detrendPolyfit = ecg - f_y;

subplot(4,1,2)
plot(tm,detrendPolyfit), grid
title(['\fontsize{16}' 'Non-linear Detrended ECG with Polyfit']), 
xlabel Sample, ylabel 'Voltage (mV)'

%% Fast Fourier Transform/Inverse-FFT
% Code from www.librow.com, by S. Chernenko
% Removes P-wave location, but can make finding R-peak easy

fresult=fft(ecg);
fresult(1 : round(length(fresult)*3/Fs)) = 0;
fresult(end - round(length(fresult)*3/Fs) : end) = 0;
detrendFFT=real(ifft(fresult));

subplot(4,1,3)
plot(tm,detrendFFT), grid
title (['\fontsize{16}' 'Detrended ECG with Fast Fourier Transform']),
xlabel Sample, ylabel 'Voltage (mV)'

%}
%% Highpass Butterworth filter
% Code from stackoverflow.com/questions/1773542/matlab-filter-noisy-ekg-signal
fNorm = 0.5/(Fs/2);    % cutoff frequency of 0.5Hz
[b,a] = butter(4, fNorm, 'high');  % 4th order filter
detrendButter = filtfilt(b, a, ecg);   % removed baseline wander

subplot(3,1,2)
plot(tm,detrendButter), grid
title (['\fontsize{16}' 'Detrended ECG with Butterworth Filter']), 
xlabel Sample, ylabel 'Voltage (mV)'

%% Powerline Interference removal
% Code from MATLAB signal processing toolbox:
%   http://uk.mathworks.com/help/signal/ug/remove-the-60-hz-hum-from-a-signal.html
if powerNoiseRemoval == true
    d = designfilt('bandstopiir','FilterOrder',2, ...
                   'HalfPowerFrequency1',58,'HalfPowerFrequency2',62, ...
                   'DesignMethod','butter','SampleRate',Fs);
    ecg = filtfilt(d,ecg);
end

subplot(3,1,3)
plot(tm,ecg), grid
title (['\fontsize{16}' 'Powerline Interference of 60Hz removed']), 
xlabel Sample, ylabel 'Voltage (mV)'






