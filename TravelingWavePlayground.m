clc; close all; clear all;
%% Load Data
load("D:\B3Raw\20231110\GangulyServer\20231110\Robot3DArrow\153829\BCI_Fixed\Data0004.mat");
fs = TrialData.Params.UpdateRate; %sampling rate
kinax = TrialData.TaskState; state1 = find(kinax==1);
state2 = find(kinax==2); state3 = find(kinax==3);
state4 = find(kinax==4);
chMap = TrialData.Params.ChMapB2;
raw = cell2mat(TrialData.BroadbandData');

%% Process Data
len_state1 = size(cell2mat(TrialData.BroadbandData(state1)'),1);
len_state2 = size(cell2mat(TrialData.BroadbandData(state2)'),1)+len_state1;
len_state3 = size(cell2mat(TrialData.BroadbandData(state3)'),1)+len_state2;

beta1 =   filter(TrialData.Params.FilterBank(4).b,...
    TrialData.Params.FilterBank(4).a,...
    raw)'; % channel*data length
beta2 = filter(TrialData.Params.FilterBank(5).b,...
    TrialData.Params.FilterBank(5).a,...
    raw)';
beta = (beta1+beta2)/2;
beta1_env =   abs(hilbert(filter(TrialData.Params.FilterBank(4).b,...
    TrialData.Params.FilterBank(4).a,...
    raw)))'; % channel*data length
beta2_env = abs(hilbert(filter(TrialData.Params.FilterBank(5).b,...
    TrialData.Params.FilterBank(5).a,...
    raw)))';
beta_env = (beta1_env+beta2_env)/2;

%% Let's see cross channel beta env
avg_beta_env = mean(beta_env);
plot(avg_beta_env)

%% Look at single channel beta burst
selected_channel = 128;
single_ch_data = beta_env(selected_channel,:);

duration = 100; % last at least 100ms
highThreshold = 0.5*(max(single_ch_data)-min(single_ch_data))+min(single_ch_data); % 80%
lowThreshold = 0.35*(max(single_ch_data)-min(single_ch_data))+min(single_ch_data); % 60%
bursts = [];
isInBurst = false;
burstStart = NaN;
for i = 1:length(single_ch_data)
    if single_ch_data(i) > highThreshold && ~isInBurst
        isInBurst = true;
        burstStart = i;
    elseif single_ch_data(i) < lowThreshold && isInBurst
        if i - burstStart >= duration
            bursts(end+1, :) = [burstStart, i];
        end
        isInBurst = false;
    end
end
if isInBurst && length(single_ch_data) - burstStart >= 100
    bursts(end+1, :) = [burstStart, length(single_ch_data)];
end

figure;
plot(raw(:,selected_channel), 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1.5,'DisplayName', 'raw channel 66');
hold on;
plot(beta(selected_channel,:), 'k', 'LineWidth', 1.5,'DisplayName', 'beta channel 66');
plot(single_ch_data, 'b','LineWidth', 1.5,'DisplayName', 'beta env channel 66');
% burst beta
for i=1:size(bursts,1)
    plot(bursts(i,1):bursts(i,2),beta(66,bursts(i,1):bursts(i,2)), 'red', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end
yline(highThreshold, '--', 'Color', 'k', 'HandleVisibility', 'off');
yline(lowThreshold, '--', 'Color', 'k', 'HandleVisibility', 'off');
legend;
line([len_state1, len_state1], ylim, 'Color', 'red', 'LineStyle', '--', 'HandleVisibility', 'off');
line([len_state2, len_state2], ylim, 'Color', 'red', 'LineStyle', '--', 'HandleVisibility', 'off');
line([len_state3, len_state3], ylim, 'Color', 'red', 'LineStyle', '--', 'HandleVisibility', 'off');


%% Look at Hg
% filter bank hg
Params=[];
Params.Fs = 1000;
Params.FilterBank(1).fpass = [70,77];   % high gamma1
Params.FilterBank(end+1).fpass = [77,85];   % high gamma2
Params.FilterBank(end+1).fpass = [85,93];   % high gamma3
Params.FilterBank(end+1).fpass = [93,102];  % high gamma4
Params.FilterBank(end+1).fpass = [102,113]; % high gamma5
Params.FilterBank(end+1).fpass = [113,124]; % high gamma6
Params.FilterBank(end+1).fpass = [124,136]; % high gamma7
Params.FilterBank(end+1).fpass = [136,150]; % high gamma8
Params.FilterBank(end+1).fpass = [30,36]; % lg1
Params.FilterBank(end+1).fpass = [36,42]; % lg2
Params.FilterBank(end+1).fpass = [42,50]; % lg3
% compute filter coefficients
for i=1:length(Params.FilterBank),
    [b,a] = butter(3,Params.FilterBank(i).fpass/(Params.Fs/2));
    Params.FilterBank(i).b = b;
    Params.FilterBank(i).a = a;
end

filtered_data=zeros(size(raw,1),size(raw,2),8);
for i=1:8 % only hg
    filtered_data(:,:,i) =  ((filter(...
        Params.FilterBank(i).b, ...
        Params.FilterBank(i).a, ...
        raw)));
end
hg_env = squeeze(mean(filtered_data.^2,3)); % hgEnvelope

%% Let's see cross channel beta env
avg_hg_env = mean(hg_env');
figure
plot(avg_hg_env)
hold on
line([len_state1, len_state1], ylim, 'Color', 'red', 'LineStyle', '--', 'HandleVisibility', 'off');
line([len_state2, len_state2], ylim, 'Color', 'red', 'LineStyle', '--', 'HandleVisibility', 'off');
line([len_state3, len_state3], ylim, 'Color', 'red', 'LineStyle', '--', 'HandleVisibility', 'off');

%% Detect burst for all channel
bursts_column = {};
for k = 1:256
    selected_channel = k;
    single_ch_data = beta_env(selected_channel,:);

    duration = 100; % last at least 100ms
    highThreshold = 0.5*(max(single_ch_data)-min(single_ch_data))+min(single_ch_data); % 80%
    lowThreshold = 0.35*(max(single_ch_data)-min(single_ch_data))+min(single_ch_data); % 60%
    bursts = [];
    isInBurst = false;
    burstStart = NaN;
    for i = 1:length(single_ch_data)
        if single_ch_data(i) > highThreshold && ~isInBurst
            isInBurst = true;
            burstStart = i;
        elseif single_ch_data(i) < lowThreshold && isInBurst
            if i - burstStart >= duration
                bursts(end+1, :) = [burstStart, i];
            end
            isInBurst = false;
        end
    end
    if isInBurst && length(single_ch_data) - burstStart >= 100
        bursts(end+1, :) = [burstStart, length(single_ch_data)];
    end
    bursts_column{k} = bursts;
end

all_index = [];
% Convert start and end to array
for i = 1:length(bursts_column)
    currentBurst = bursts_column{i};
    for k = 1:size(currentBurst,1)
        index_list = currentBurst(k,1):currentBurst(k,2);
        all_index = [all_index, index_list];
    end
end

index_count = accumarray(all_index', 1);
index_frequency = index_count / 253;
figure;
bar(1:numel(index_count), index_count);
xlabel('Index');
ylabel('Frequency');
title('Distribution of Indices');

%% Cross exam avg_beta_env and beta_burst_count
figure
yyaxis left;
plot(avg_beta_env)
ylabel('average beta envolope')
hold on
yyaxis right;
plot(index_frequency)
ylabel('Beta burst frequency')

% Look at the correlation coefficient
beta_env_burst_corre = corrcoef(avg_beta_env, index_frequency);

%% relocate beta structrue to 3D
% Initialize the new data matrix
newBetaEnv = NaN(11, 23, length(beta));  % Use NaN to easily identify unmapped elements

% Iterate through each element of chMap to re-coordinate the data
for row = 1:size(chMap, 1)
    for col = 1:size(chMap, 2)
        channelNumber = chMap(row, col);  % Get the channel number from chMap
        if ~isnan(channelNumber) && channelNumber <= 256 && channelNumber >= 1  % Check if the channelNumber is valid
            % Copy the data from this channel to the new location
            newBetaEnv(row, col, :) = beta(channelNumber, :);
        end
    end
end

%% animate it
% Define the figure
figure;

% Number of time points
numTimePoints = size(newBetaEnv, 3);

% Loop through each time point to update the plot
for t = 2000:numTimePoints
    % Update the image data with the slice at the current time point
    imagesc(squeeze(newBetaEnv(:,:,t)));
    
    % Optional: Set the color axis if you want consistent coloring across frames
    caxis([min(newBetaEnv,[],'all') max(newBetaEnv,[],'all')]);
    
    % Add a title with the current time point (optional)
    title(sprintf('Time: %d', t));
    
    % Add colorbar and adjust its limits to match your data range (optional)
    colorbar;
    
    % Pause briefly to create the animation effect
    pause(0.01);  % Adjust the pause time as needed for smooth playback
end

%% IDing of oscillation using Morlet wavelets

raw = raw';
% Define frequency range
frequencies = logspace(log10(3), log10(40), 200); % 200 frequencies from 3 to 40 Hz

% Number of cycles in the wavelet, typically 3-10, affects time/frequency resolution
nCycles = 6;

% Pre-allocate power matrix
power = zeros(length(frequencies), size(raw, 1), size(raw, 2));

% Loop over channels
for channel = 1:size(raw, 1)
    signal = raw(channel, :); % Signal for the current channel

    % Loop over frequencies
    for fIdx = 1:length(frequencies)
        f = frequencies(fIdx); % Current frequency

        % Generate Morlet wavelet
        st = 1 / (2 * pi * (f / nCycles));
        t = -3.5 * st : 1/fs : 3.5 * st;
        wavelet = exp(2 * 1i * pi * f .* t) .* exp(-t.^2 / (2 * st^2));

        % Convolve signal with wavelet
        signalFFT = fft(signal, length(signal) + length(wavelet) - 1);
        waveletFFT = fft(wavelet, length(signal) + length(wavelet) - 1);
        convolutionResult = ifft(signalFFT .* waveletFFT);

        % Trim the convolutionResult to match the original signal length more accurately
        convolutionLength = length(convolutionResult);
        expectedLength = length(signal);
        startIdx = floor((convolutionLength - expectedLength) / 2) + 1;
        endIdx = startIdx + expectedLength - 1;

        % Adjust if the endIdx exceeds the convolutionLength (rare edge case)
        if endIdx > convolutionLength
            endIdx = convolutionLength;
            startIdx = endIdx - expectedLength + 1;
        end

        convolutionResult = convolutionResult(startIdx:endIdx);

        % Confirm the convolutionResult is the correct length
        assert(length(convolutionResult) == length(signal), ...
               ['Mismatch in convolution result length for channel %d, ' ...
                'frequency %.2f Hz. Expected %d, got %d'], ...
                channel, f, length(signal), length(convolutionResult));

        % Compute power and store
        power(fIdx, channel, :) = abs(convolutionResult).^2;
    end
end

%%
meanPowerSpectrum = mean(power, 3); % Mean across time (3rd dimension of power matrix)
meanPowerSpectrum = squeeze(meanPowerSpectrum);

% Fit a line in Log-Log Coordinates
logFreq = log10(frequencies); % Log-transform frequencies
coefficients = cell(size(meanPowerSpectrum, 1), 1); % Pre-allocate for coefficients
for channel = 1:size(meanPowerSpectrum, 2)
    logPower = log10(meanPowerSpectrum(:, channel)); % Log-transform power for this channel
    
    % Perform robust linear regression
    [b, ~] = robustfit(logFreq, logPower);
    
    % Calculate the fit line and subtract it from the actual power spectrum
    fitLine = b(1) + b(2) * logFreq;
    normalizedPowerSpectrum(:, channel) = logPower - fitLine';
end

%%
peaksInfo = cell(1, size(meanPowerSpectrum, 2)); % Pre-allocate for peaks information

for channel = 1:size(meanPowerSpectrum, 2)
    [pks, locs] = findpeaks(normalizedPowerSpectrum(:, channel));
    
    % Calculate mean and standard deviation for significant peaks detection
    meanPeak = mean(normalizedPowerSpectrum(:, channel));
    stdPeak = std(normalizedPowerSpectrum(:, channel));
    
    % Identify significant peaks
    significantPeaksIdx = pks > (meanPeak + stdPeak); % Peaks more than 1 std above mean
    significantPeaks = pks(significantPeaksIdx);
    significantFreqs = frequencies(locs(significantPeaksIdx));
    
    % Store significant peaks information
    peaksInfo{channel} = struct('Peaks', significantPeaks, 'Frequencies', significantFreqs);
end

%% Plot them
% Define the channel of interest
channelOfInterest = 1;

% Select the data for the channel of interest
meanPowerChannel = meanPowerSpectrum(:, channelOfInterest);
logMeanPowerChannel = log10(meanPowerChannel);

% Perform the robust linear regression on log-log data
[b, stats] = robustfit(log10(frequencies), logMeanPowerChannel);
fitLine = b(1) + b(2) * log10(frequencies);

% Plot figure D: Log Power Spectrum with Regression Line
figure;
loglog(frequencies, meanPowerChannel, 'b'); % Log power in blue
hold on;
loglog(frequencies, 10.^fitLine, 'k--'); % Regression line in black dashed
xlabel('Frequency (Hz)');
ylabel('Log Power (µV^2/Hz)');
title('Log Power Spectrum with 1/f Trend');
hold off;

% Calculate z-scores for the normalized power
zScores = (normalizedPowerSpectrum(:, channelOfInterest) - mean(normalizedPowerSpectrum(:, channelOfInterest))) ./ std(normalizedPowerSpectrum(:, channelOfInterest));

% Plot figure E: Normalized Power Spectrum as Z-Score
figure;
semilogx(frequencies, zScores, 'r'); % Z-scores in red
xlabel('Frequency (Hz)');
ylabel('Z-scored Power');
title('Normalized Power Spectrum as Z-Score');
% Horizontal line at Z = 1
hold on;
line(get(gca,'xlim'), [1 1], 'Color', 'black', 'LineStyle', '--');
hold off;


%% Sanity check by only using PSD
[pxx,frequencies] = pwelch(raw(1,:), [], [], [], 1000);

figure
plot(frequencies, 10*log10(pxx))
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
title('Power Spectral Density of Signal')

% Define frequency range
freqRange = (frequencies >= 3) & (frequencies <= 40);

figure
% Plot only the selected frequency range
plot(log10(frequencies(freqRange)), 10*log10(pxx(freqRange)))
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
title('Power Spectral Density of Signal (3-40 Hz)')

%% Look at four connecting ECoG Channel and relationship
ch_num = [114, 109, 76, 72];
% Assuming beta is your data matrix
beta_ch_num = beta(ch_num,:); % Your original line
beta_ch_num_look = beta_ch_num(:,10:1010); % Select a subset of your data

x = 1:size(beta_ch_num_look,2); % X-axis values

figure; % Open a new figure window
hold on; % Keep the plot window active to overlay multiple plots

% Define colors for each signal line
colors = ['r', 'g', 'b', 'm']; % Example colors: red, green, blue, magenta

% Define offsets manually or compute them automatically
offsets = [0, 1, 2, 3]; % Example offsets, adjust as needed

% Store handles for the main data lines for the legend
h = zeros(1, length(ch_num)); % Initialize array for handles

% Plot each row with an offset and mark peaks with vertical dashed lines
for i = 1:size(beta_ch_num_look,1)
    y_data = beta_ch_num_look(i,:) + offsets(i); % Adjusted Y data with offset
    h(i) = plot(x, y_data, 'LineWidth', 2, 'Color', colors(i)); % Plot the data with specified color and store handle
    
    % Find peaks in the adjusted data
    [peaks, locations] = findpeaks(y_data);
    
    % Get current y-axis limits
    yLimits = ylim;
    
    % Mark the peaks with vertical dashed lines across the entire y-axis range
    for j = 1:length(peaks)
        line([x(locations(j)) x(locations(j))], [yLimits(1), yLimits(2)], 'LineStyle', '--', 'Color', colors(i));
    end
end

% Use the handles to create a legend for the main lines only
legend(h, 'ch 114', 'ch 109', 'ch 76', 'ch 72');
hold off; % No more plots to add

%% look at gradient
% Assuming your data is in a 3D matrix called 'data' of size 11x23xn

% Get the size of your data
[rows, cols, nTimeSteps] = size(newBetaEnv);

% Preallocate arrays to store gradients (for efficiency)
gradX = zeros(rows, cols, nTimeSteps);
gradY = zeros(rows, cols, nTimeSteps);

for t = 1:nTimeSteps
    [GX, GY] = gradient(newBetaEnv(:, :, t));
    gradX(:, :, t) = GX;
    gradY(:, :, t) = GY;
    
    % Optional: Visualize the gradient magnitude for this time step
    gradMag = sqrt(GX.^2 + GY.^2); % Calculate the magnitude of the gradient

end

% Prepare for animation
figure;
[X, Y] = meshgrid(1:cols, 1:rows);

% Initialize the animation with the first set of vectors
quiverHandle = quiver(X, Y, gradX(:,:,1), gradY(:,:,1), 'AutoScale', 'off');
axis tight; % Fix the axes limits
title('Gradient Vectors - Time Step 1');

% Setup VideoWriter
videoFileName = 'gradientVectorsAnimation.avi';
v = VideoWriter(videoFileName);
open(v);

% Animate over time steps and save to video
for t = 1:nTimeSteps
    set(quiverHandle, 'UData', gradX(:,:,t), 'VData', gradY(:,:,t));
    title(['Gradient Vectors - Time Step ', num2str(t)]);
    drawnow; % Update the plot
    
    % Capture the frame
    frame = getframe(gcf); % If you're using a specific figure or axes, use getframe(figureHandle) or getframe(axesHandle)
    writeVideo(v, frame); % Write the frame to the video
end

% Clean up
close(v); % Close the VideoWriter object

%% Look at phase syncrhonization analysis
% Example for a single time series
signal114 = newBetaEnv(9,8,:); % Extract a single time series for demonstration
analyticSignal114 = hilbert(signal114);
instantaneousPhase114 = angle(analyticSignal114); % Phase of the analytic signal
figure
plot(squeeze(instantaneousPhase114));

signal109 = newBetaEnv(9,9,:); % Extract a single time series for demonstration
analyticSignal109 = hilbert(signal109);
instantaneousPhase109 = angle(analyticSignal109); % Phase of the analytic signal
figure
plot(squeeze(instantaneousPhase109));

% Assuming instantaneousPhase1 and instantaneousPhase2 are the phases of two signals
phaseDifference = squeeze(instantaneousPhase114) - squeeze(instantaneousPhase109);
plv = abs(mean(exp(1i * phaseDifference))); % Compute PLV

%% Look at phase velocity analysis

% Initialize the matrix to store the instantaneous phase
instantaneousPhase = zeros(size(newBetaEnv)); % data is your original 11x23xn dataset

for i = 1:11
    for j = 1:23
        signal = squeeze(newBetaEnv(i, j, :)); % Extract the time series for each spatial point
        analyticSignal = hilbert(signal);
        instantaneousPhase(i, j, :) = angle(analyticSignal); % Compute the instantaneous phase
    end
end

% Preallocate arrays for the phase velocity components
phaseVelocityX = zeros(size(newBetaEnv));
phaseVelocityY = zeros(size(newBetaEnv));

for t = 1:4800
 % Assuming 'n' is the number of time steps
    [GX, GY] = gradient(squeeze(instantaneousPhase(:, :, t)));
    phaseVelocityX(:, :, t) = GX;
    phaseVelocityY(:, :, t) = GY;
end

% Calculate magnitude and direction (example for a single time step)
t = 1; % Example time step
magnitude = sqrt(phaseVelocityX(:,:,t).^2 + phaseVelocityY(:,:,t).^2);
direction = atan2(phaseVelocityY(:,:,t), phaseVelocityX(:,:,t));

% Visualization (example for the same single time step)
[X, Y] = meshgrid(1:23, 1:11); % Assuming the second dimension is 'x' and the first is 'y'
% Initialize the figure for animation
fig = figure;
set(fig, 'Position', [100, 100, 600, 450]); % Set the figure window size [left bottom width height]

% Create the initial quiver plot
axesHandle = axes('Parent', fig, 'NextPlot', 'replaceChildren');
quiverHandle = quiver(axesHandle, X, Y, phaseVelocityX(:,:,1), phaseVelocityY(:,:,1), 'AutoScale', 'off');

% Lock the axes limits and remove ticks and labels
xlim(axesHandle, [1 cols]);
ylim(axesHandle, [1 rows]);
set(axesHandle, 'XTick', [], 'YTick', [], 'XTickLabel', {}, 'YTickLabel', {});
axis(axesHandle, 'manual'); % Lock the axis limits
axis(axesHandle, 'off'); % Turn off the axis

% Keep the aspect ratio consistent
axis(axesHandle, 'equal');

title('Phase Velocity Vectors - Time Step 1');
xlabel('X-axis');
ylabel('Y-axis');

% Prepare for video writing
videoFileName = 'PhaseVelocityAnimation.avi';
v = VideoWriter(videoFileName);
open(v);

% Define a maximum vector length
maxLength = 2; % Example value, adjust this to fit your data appropriately

% Animate through the time steps and cap the vector lengths
for t = 1:nTimeSteps
    U = phaseVelocityX(:,:,t);
    V = phaseVelocityY(:,:,t);
    vectorLengths = sqrt(U.^2 + V.^2);

    % Find vectors longer than the maximum length
    longVectors = vectorLengths > maxLength;
    
    % Cap the vectors at the maximum length while preserving the direction
    U(longVectors) = U(longVectors) ./ vectorLengths(longVectors) * maxLength;
    V(longVectors) = V(longVectors) ./ vectorLengths(longVectors) * maxLength;
    
    % Update the quiver plot
    set(quiverHandle, 'UData', U, 'VData', V);
    title(sprintf('Phase Velocity Vectors - Time Step %d', t));
    drawnow;
    
    % Capture the frame for video
    frame = getframe(fig); % Capture the figure window
    writeVideo(v, frame);
end

% Close the video file
close(v);

% Notify user
disp(['Animation saved to ', videoFileName]);

