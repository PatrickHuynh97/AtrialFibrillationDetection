%% Version 2
% Operates on a defined width, fixed window of a signal
% Iterative operation, simulating real life function where data is slowly captured.
clear all
addpath('mcode');
addpath('database');

%% Define parameters for optimal performance
mV = -10; 
rWaveCutoff = 0.45;
p_wave_window_start = 0.2; % how far to look behind R-wave peak for P-wave

lead = 1; % which ECG lead to analyse 

butterworth = true; % baseline wander removal
powerNoiseRemoval = true; % Powerline Interference removal (60hz)

windowSize = 16; % Optimal window size
signalStart = 1; % Default signal starting point
signalTime = 100; % Amount of signal to download (in seconds)

sDev_threshold = 0.2; % Optimal SDSD threshold
p_score_threshold = 10; % Optimal Pscore threshold

plotAllWindows = false; % Plot all windows
plotAFWindows = false; % Plot only windows exceeding AF thresholds

%% Training data:

% AF
signalCode = 'afdb/04015'; % requires butterworth

%{
signalCode = 'afdb/04048'; % requires butterworth and powerline 
signalTime = 40;
%}


% NSR 
%{
signalCode = 'nsrdb/16265'; 
[s,Fs,t]=rdsamp(signalCode,[], 250*1); % Load 1 second to get sample rate
signalStart = 1*250;
%}

%{
signalCode = 'mitdb/103'; % NSR at 1:09. 
%}

% Other Arrhythmia
%{
signalCode = 'mitdb/106'; % requires butterworth
%}

%{
signalCode = 'svdb/805'; 
%}

%% Testing data

%{
signalCode = 'afdb/05091';
%}

%% Downloading signal 
% Get samplerate, download signalTime seconds of signal

[s,Fs,t]=rdsamp(signalCode,[], 250*1); % Load 1 second to get sample rate

signalEnd = signalStart*Fs + (Fs*signalTime); % Locate end of signal 
[signal,Fs,tm]=rdsamp(signalCode,[], signalEnd, signalStart); 
ecg = signal(:,lead); % select lead defined in parameters

%% Optional Butterworth High Pass Filter to remove baseline wander and Powerline Interference removal
% 	Code: stackoverflow.com/questions/1773542/matlab-filter-noisy-ekg-signal

if butterworth == true
    fNorm = 1.1/(Fs/2);                 % cutoff frequency of 0.5Hz
    [b,a] = butter(6, fNorm, 'high');   % 6th order filter
    ecg = filtfilt(b, a, ecg);
end

% Code from MATLAB signal processing toolbox:
% 	http://uk.mathworks.com/help/signal/ug/remove-the-60-hz-hum-from-a-signal.html
if powerNoiseRemoval == true
    d = designfilt('bandstopiir','FilterOrder',2, ...
                   'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
                   'DesignMethod','butter','SampleRate',Fs);
    ecg = filtfilt(d,ecg);
end


%% Detecting R-wave peaks

r_wave_peak = zeros(1,10000);      % r-wave peak position
bpw = 0; % number of beats detected

% R-wave peak detection
for k=1:length(ecg)         
    % If point has an amplitude > threshold and is a turning point we should consider it 
    if(ecg(k) > rWaveCutoff && ecg(k) > ecg(k-2) && ecg(k) > ecg(k+1))
        % Obtain average of the gradient between k following 3 samples
        gradients = zeros(1,3);
        for i=1:3
            yDif = ecg(k+i) -  ecg(k+i+1);  
            xDif = tm(k+i) - tm(k+i+1);   
            gradients(i) = yDif / xDif;  
        end
        
        meanGradient = mean(gradients);
        
        if (meanGradient < mV)  % gradient should be negative and steep according to defined parameters
            % make sure we only check for double detections after first run
            if r_wave_peak(1) ~= 0
                closeness = tm(k) - tm(r_wave_peak(bpw)); 
                if closeness > 0.1 % if the detected peaks are further than 0.1s apart
                    bpw = bpw + 1; % store for easy BPM calculation
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

%% Detecting P-wave peaks
p_wave_peak = zeros(1,1000);       % p-wave peak position
p_wave_window = zeros(1000,2);     % start/end of p-wave 
p_halftime = round(Fs*0.06, 0);     % half the duration of an average p-wave

% window in which P-wave must occur (200ms before r-wave peak)
for k=1:length(r_wave_peak)
    p_wave_window(k,:) = [r_wave_peak(k) - round(p_wave_window_start*Fs, 0) r_wave_peak(k) - round(Fs*0.05, 0)];
end

p_wave_window(all(~p_wave_window,2), : ) = []; % remove zero values
p_wave_window = p_wave_window(all(p_wave_window > 0, 2), : ); % ignore negative rows

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
    
    % set new p-wave window according to 0.12s normal duration and new peak
    p_wave_window(k,:) = [(p_wave_peak_tm - p_halftime), (p_wave_peak_tm + p_halftime)];
end

p_wave_peak(p_wave_peak == 0) = [];

%% Simulate running in real time 
% Establish moving window of size windowSize
% Calculate SDev and P-score, store within each window

totalIterations = signalTime - windowSize; % Total number of windows
AF_score = zeros(totalIterations,2); % Store SDev and P-Score
plots = 0; % control how many windows are plotted. 0 = all windows are plotted

p_wave_score = zeros(1,length(p_wave_peak));   % store all scores in matrix

% Calculate Standard Deviation and P-Score within fixed window
for i=1:(totalIterations+1)   
    %% RR-interval standard deviation    
    % On first iteration, start with very first sample
    if i==1
        window_start = 1;
    else 
        window_start = (i-1) * Fs;
    end
    window_end = window_start + (Fs * windowSize);
    
    % find R peaks in the window we are considering
    consider = find(r_wave_peak > window_start & r_wave_peak < window_end);
    r_wave_differences = zeros(1,length(consider-1));
    
    % Calculate Difference in consecutive R-waves
    for k=1: (length(consider)-1)
        r_wave_differences(k) = (tm(r_wave_peak(consider(k+1))) - tm(r_wave_peak(consider(k))))^2;
    end

    standard_deviation = std(r_wave_differences); 
    AF_score(i,1) = standard_deviation; 
    
    
    %% P-wave Score
    
    for l=1:length(consider) 
        
        p_wave_peak_dif_before = ecg(p_wave_peak(consider(l))) - ecg(p_wave_window(consider(l),1));
        p_wave_peak_dif_after = ecg(p_wave_peak(consider(l))) - ecg(p_wave_window(consider(l),2));
        
        p_score = abs(1/p_wave_peak_dif_before - 1/p_wave_peak_dif_after);

        % if p_score is inf, set it to a very small value instead
        if p_score == inf 
            p_score = 0.00001;
        end
        
        p_wave_score(l) =  p_score;
        
        p_wave_score(p_wave_score == 0) = []; 
        mean_p_wave = mean(p_wave_score);
        AF_score(i,2) = mean_p_wave;
    end
    %% Plots
    % plot every window where AF is diagnosed 
    if plotAFWindows == true 
        if mean_p_wave > p_score_threshold && standard_deviation > sDev_threshold
           plots = plots +1;
           if plots < 20
                hold on
                    plots = plots + 1;
                    figure()
                    plot(tm(window_start : window_end), ecg(window_start : window_end), ...
                        tm(r_wave_peak(consider)), ecg(r_wave_peak(consider)), 'ro', ...
                        tm(p_wave_peak(consider)), ecg(p_wave_peak(consider)), 'go')
                    xlabel 'Sample (s)', ylabel 'Voltage (mV)'
                    titleString = ['Sliding window from ', num2str(window_start/Fs), ' seconds to ', ...
                        num2str(window_end/Fs),' seconds where SDev = ', num2str(standard_deviation), ...
                        ' and mean P-score = ', num2str(mean_p_wave)];
                    title (['\fontsize{12}' titleString])
                hold off
           end
        end
        
    elseif plotAllWindows == true
        
        hold on
            figure()
            plot(tm(window_start : window_end), ecg(window_start : window_end), ...
                tm(r_wave_peak(consider)), ecg(r_wave_peak(consider)), 'ro', ...
                tm(p_wave_peak(consider)), ecg(p_wave_peak(consider)), 'go')
            xlabel 'Sample (s)', ylabel 'Voltage (mV)'
            titleString = ['Sliding window from ', num2str(window_start/Fs), ' seconds to ', ...
                num2str(window_end/Fs),' seconds where SDev = ', num2str(standard_deviation), ...
                ' and mean P-score = ', num2str(mean_p_wave)];
            title (['\fontsize{12}' titleString])
        hold off
      
    elseif i == 1 % if both plotting parameters are false, just plot the entire signal
        
        hold on
            plot(tm, ecg, ...
                tm(r_wave_peak), ecg(r_wave_peak), 'ro', ...
                tm(p_wave_peak), ecg(p_wave_peak), 'go');
      
            xlabel 'Sample (s)', ylabel 'Voltage (mV)'
         
            titleString = [num2str(signalTime), ' seconds of ', signalCode];
            title (['\fontsize{12}' titleString])
        hold off
        
    end 
end
