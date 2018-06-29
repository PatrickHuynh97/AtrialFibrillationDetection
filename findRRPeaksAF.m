%% Uses MIT-BIH AF Database and identifies RR peaks and RR-Interval (heart-rate)
%   1. Applies a Butterworth High Pass filter to remove Baseline Wander.
%   2. Apply R-wave detection algorithm
clear all
addpath('mcode');

mV = -20;
rWaveCutoff = 0.45;

% Download ECG Signal from Atrial Fib. Database, select lead 1
[signal,Fs,tm]=rdsamp('afdb/04015',[], 250*100); 
ecgNoisy = signal(:,1);    

% Testing with Normal Sinus Rhythm 
%[signal,Fs,tm]=rdsamp('nsrdb/16265',[],12800); 
%ecg = signal(:,1);    

%% Butterworth High Pass Filter to remove baseline wander

fNorm = 0.5/(Fs/2);                 % cutoff frequency of 0.5Hz
[b,a] = butter(6, fNorm, 'high');   % 6th order filter
ecg = filtfilt(b, a, ecgNoisy);

%% Locate R-wave peaks

bpw = 0;                    % keeps running total of peaks located in window
r_wave_ecg = zeros(1,30);   % Prefill array for efficiency
r_wave_tm = zeros(1,30);    % Prefill array for efficiency

for k=1:length(ecg)         % loop through all samples
    % If point has an amplitude > 0.45mV  and is a turning point we should consider it 
    if(ecg(k) > rWaveCutoff && ecg(k) > ecg(k-1) && ecg(k) > ecg(k+1))
              
        % Obtain average of the gradient between k following 3 samples
        gradients = zeros(1,3);
        for i=1:3
            yDif = ecg(k+i) -  ecg(k+i+1);  
            xDif = tm(k+i) - tm(k+i+1);   
            gradients(i) = yDif / xDif;  
        end
        
        meanGradient = mean(gradients);
        
        if (meanGradient < mV)  % gradient should be negative and steep 
            bpw = bpw + 1;              % store for easy BPM calculation
            r_wave_ecg(bpw) = ecg(k) ;  % store amplitude of peak 
            r_wave_tm(bpw) = tm(k);     % store time of peak
        end
    end
end

% Remove preallocated 0's from arrays, calculate BPM
r_wave_ecg(r_wave_ecg == 0) = [];
r_wave_tm(r_wave_tm == 0) = [];
bpm = bpw * (60/18);

% Plot isolated R-waves and corresponding peaks
hold all
subplot(1,1,1)
plot(tm, ecg,r_wave_tm, r_wave_ecg, 'ro'),
xlabel 'Sample (s)', ylabel 'Voltage (mV)'
title 'R-wave Peak Detection'













