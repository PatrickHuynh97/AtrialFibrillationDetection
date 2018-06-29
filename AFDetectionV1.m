%% Version 1
% Operates on an entire signal
clear all
addpath('mcode');

%% Downloading Signal 
mV = -20;
rWaveCutoff = 0.45;

butterworth = true;

% Download ECG Signal from Atrial Fib. Database, select lead 1
[signal,Fs,tm]=rdsamp('afdb/04015',[], 250*100); % 08215 is also good
ecg = signal(:,1);    

%Testing with Normal Sinus Rhythm (comment out butterworth)
%[signal,Fs,tm]=rdsamp('nsrdb/16265',[],128*500); 
%ecg = signal(:,1);    

%% Butterworth High Pass Filter to remove baseline wander
if butterworth == true
    fNorm = 0.5/(Fs/2);                 % cutoff frequency of 0.5Hz
    [b,a] = butter(6, fNorm, 'high');   % 6th order filter
    ecg = filtfilt(b, a, ecg);
end
%% R-wave 

r_wave_peak = zeros(1,1000);      % r-wave peak position
bpw = 0; % number of beats detected

% R-wave peak detection
for k=1:length(ecg)         
    % If point has an amplitude > 0.45mV  and is a turning point we should consider it 
    if(ecg(k) > rWaveCutoff && ecg(k) > ecg(k-2) && ecg(k) > ecg(k+1))
              
        % Obtain average of the gradient between k following 3 samples
        gradients = zeros(1,3);
        for i=1:3
            yDif = ecg(k+i) -  ecg(k+i+1);  
            xDif = tm(k+i) - tm(k+i+1);   
            gradients(i) = yDif / xDif;  
        end
        
        meanGradient = mean(gradients);
        
        if (meanGradient < mV)  % gradient should be negative and steep 
            % make sure we only check for double detections after first run
            if r_wave_peak(1) ~= 0
                closeness = tm(k) - tm(r_wave_peak(bpw)); 
                if closeness > 0.1 % if the detected peaks are further than 0.1s apart
                    bpw = bpw + 1;              % store for easy BPM calculation
                    r_wave_peak(bpw) = k ;  % store array position of peak
                end
                % otherwise, ignore the detection
            else % first time code runs, take first R-wave peak found as correct
                bpw = bpw + 1;              % store for easy BPM calculation
                r_wave_peak(bpw) = k ;  % store array position of peak
            end
        end
    end
end
 
r_wave_peak(r_wave_peak == 0) = []; % remove zero values

%% RMSSD and SDNN calculation
length(r_wave_peak)
r_wave_differences = zeros(1,length(r_wave_peak)-1);

% loop through all peaks, calculate differences between one and the next
for i=1: (length(r_wave_peak)-1)
    r_wave_differences(i) = (tm(r_wave_peak(i+1)) - tm(r_wave_peak(i))); 
end

rmssd = sqrt(1/length(r_wave_differences) * sum(r_wave_differences.^2));

sdev = std(r_wave_differences);

%% P-wave 
p_wave_peak = zeros(1,10000);      % p-wave peak position
p_wave_window = zeros(10000,2);    % start/end of p-wave 
p_wave_gradients = zeros(0,2);  % before/after gradients of peak  
p_halftime = round(Fs*0.06, 0); % distance between middle and end of p-wave

% window in which P-wave must occur (200ms before r-wave peak)
for k=1:length(r_wave_peak)
    p_wave_window(k,:) = [r_wave_peak(k) - round(0.2*Fs, 0) r_wave_peak(k) - round(Fs*0.05, 0)];
end

p_wave_window( all(~p_wave_window,2), : ) = []; % remove zero values

% Search for peak point in each P-wave
for k=1:length(p_wave_window)
    
    p_wave_peak_value = -10; % initialise value to store peak
    p_wave_peak_tm = 0;
    
    % Loop through P-wave window
    for p = p_wave_window(k,1): p_wave_window(k,2)
        if ecg(p) > p_wave_peak_value % if the current sample is greater than current peak
            p_wave_peak_value = ecg(p);    % save the new peak
            p_wave_peak_tm = p;
        end
    end
    p_wave_peak(k) = p_wave_peak_tm;
    
    % set new p-wave window according to 0.12s normal duration
    p_wave_window(k,:) = [(p_wave_peak_tm - p_halftime), (p_wave_peak_tm + p_halftime)];
end

p_wave_peak(p_wave_peak == 0) = [];

% loop through all windows, calculate gradients
for k=1:length(p_wave_peak)
    y_dif_backward = ecg(p_wave_peak(k)) - ecg(p_wave_window(k, 1));
    x_dif_backward = tm(p_wave_peak(k)) - tm(p_wave_window(k, 1));
    
    y_dif_forward = ecg(p_wave_window(k, 2)) - ecg(p_wave_peak(k));
    x_dif_forward = tm(p_wave_window(k, 2)) - tm(p_wave_peak(k));
    
    p_wave_gradients(k,:) = [y_dif_backward/x_dif_backward y_dif_forward/x_dif_forward];
end

p_wave_gradients( all(~p_wave_gradients,2), : ) = [];

%% Implementing Version 2
% The initial and final amplitudes should be closer to each other than peak
% Strength of P-wave will be stored in p_wave_scores
p_wave_amount = length(p_wave_gradients); % number of P-waves detected
p_wave_amp_dif = zeros(1,p_wave_amount);
p_score_v2 = zeros(1,p_wave_amount);
p_wave_gradient_dif = zeros(1,p_wave_amount);

graph_size = 100;

% Calculate difference (absolute) in amplitude between peaks and start/end
for i=1:p_wave_amount
    p_wave_amp_dif(i) = abs(ecg(p_wave_window(i,2)) + ecg(p_wave_window(i,1))); 
    p_wave_gradient_dif(i) = abs(p_wave_gradients(i,1) + p_wave_gradients(i,2));
    p_wave_peak_dif_before = ecg(p_wave_peak(i)) - ecg(p_wave_window(i,1));
    p_wave_peak_dif_after = ecg(p_wave_peak(i)) - ecg(p_wave_window(i,2));
    
    p_score_v2(i) = abs((1/p_wave_peak_dif_before) - (1/p_wave_peak_dif_after));
    
    if p_score_v2(i) > 10 && (p_wave_peak(i)-graph_size) > 0 && (p_wave_peak(i) + graph_size) <length(tm)
        % Plot the section of the signal which passed the threshold  
        %{
        figure()
        hold on
        plot(tm((p_wave_peak(i)-graph_size):(p_wave_peak(i)+graph_size)), ...
        ecg((p_wave_peak(i)-graph_size):(p_wave_peak(i)+graph_size)))
        plot(tm(r_wave_peak(i)), ecg(r_wave_peak(i)), 'ro')
        plot(tm(p_wave_peak(i)), ecg(p_wave_peak(i)), 'go')
        hold off
        %}
    end
end

%% Implementing Version 3 / Parson's P-wave Detection
% Calculate average amplitude in each P-wave window
% Sum the following values:
%   Absolute difference between average and peak 
%   Absolute difference between initial/final amp. and average
% The higher the value, the stronger the P-wave

p_wave_amount = length(p_wave_gradients); % number of P-waves detected
p_score_v3 = zeros(1,p_wave_amount); % p_score for each wave
graph_size = 300;

for i=1:p_wave_amount
    
    % Calculate average amplitude in P-wave window
    average = mean(ecg(p_wave_window(i,1):p_wave_window(i,2)));
    
    initial_to_average = abs(average - ecg(p_wave_window(i,1)));
    end_to_average = abs(average - ecg(p_wave_window(i,2)));
    peak_to_average = abs(average - ecg(p_wave_peak(i)));
    
    p_score_v3(i) = initial_to_average + end_to_average + peak_to_average;
    
    if p_score_v3(i) > 0.5 && (p_wave_peak(i)-graph_size) > 0 && (p_wave_peak(i) + graph_size) <length(tm)
        % Plot the section of the signal which passed the threshold
        
        figure()
        hold on
        plot(tm((p_wave_peak(i)-graph_size):(p_wave_peak(i)+graph_size)), ...
        ecg((p_wave_peak(i)-graph_size):(p_wave_peak(i)+graph_size)))
        plot(tm(r_wave_peak(i)), ecg(r_wave_peak(i)), 'ro')
        plot(tm(p_wave_peak(i)), ecg(p_wave_peak(i)), 'go')
        hold off
        
    end
end


%% Plot
plot(tm, ecg,tm(r_wave_peak), ecg(r_wave_peak), 'ro', tm(p_wave_peak), ecg(p_wave_peak), 'go')
for k=1:length(p_wave_window)
    line([tm(p_wave_window(k, 1)) tm(p_wave_window(k, 1))],[-1 3], 'Color', 'black')
    line([tm(p_wave_window(k, 2)) tm(p_wave_window(k, 2))],[-1 3], 'Color', 'black')
end
xlabel 'Sample (s)', ylabel 'Voltage (mV)'
title 'P-wave analysis'
