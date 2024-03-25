close all; clear all; clc
%%
load("D:\B3Raw\20231110\GangulyServer\20231110\Robot3DArrow\153829\BCI_Fixed\Data0004.mat");
% file = 'G:\Ganguly lab\data\B1\20220902\GangulyServer\20220902\Robot3DArrow\100727\BCI_Fixed\Data0008.mat';
% load(file)
raw = cell2mat(TrialData.BroadbandData');
chMap = TrialData.Params.ChMap;
Fs = 1000; % 1000 Hz
bpFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',3,'HalfPowerFrequency2',40, ...
    'SampleRate',1e3);
peakFreqs = zeros(1, 128);
for i=1
% data = filtfilt(bpFilt,raw(:,i));
data = raw(:,i);
[wt,f] = cwt(data, 'amor', Fs); % 'amor' - Morlet
power = mean((abs(wt).^2),2);
logf = log10(f);
logPower = log10(power);
% fit
% [bhat p wh se ci t_stat]=robust_fit(logf,logPower,1);
figure
plot(logf,logPower)
% hold on; plot(logf,bhat(1)+bhat(2).*logf,'r')
% figure;
% plot(logf,logPower - bhat(1)+bhat(2).*logf)
    % Detrended log power
%     detrendedLogPower = logPower - (bhat(1) + bhat(2).*logf);
    % Find peaks
%     [pks, locs] = findpeaks(detrendedLogPower, 'MinPeakProminence', 0.1);
%     % Check if peaks were found
%     if ~isempty(pks)
%         % Find the frequency of the largest peak
%         [~, maxIdx] = max(pks);
%         peakFreqs(i) = round(f(locs(maxIdx)));
%     else
%         peakFreqs(i) = NaN;
%     end
end