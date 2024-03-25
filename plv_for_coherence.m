% % main
% 
% 
% clc
% clear all
% close all
% 
% clc; close all; clear all
% 
% 
% root_path='G:\Ganguly lab\data\B1\';
% % foldernames = {'20220608', '20220622', '20220624', '20220629','20220701', '20220713', '20220715','20220720', '20220722','20220727', '20220803', '20220810', '20220817','20220831', '20220902'};
% 
% foldernames = {'20210613','20210616','20210623','20210625','20210630','20210702',...
%     '20210707','20210716','20210728','20210804','20210806','20210813','20210818',...
%     '20210825','20210827','20210901','20210903','20210910','20210917','20210924','20210929',...
%     '20211001''20211006','20211008','20211013','20211015','20211022','20211027','20211029','20211103',...
%     '20211105','20211117','20211119','20220126','20220128','20220202','20220204','20220209','20220211',...
%     '20220218','20220223','20220225','20220302','20220309','20220311',...
%     '20220316','20220323','20220325','20220520','20220722','20220727','20220729',...
%     '20220803','20220810','20220812','20220817','20220819','20220831','20220902',...
%     '20220907'};
% 
% cd(root_path)
% 
% task_files={}; python_files = [];ph_all = []; phi_all = []; r1_all =[]; vel_all = [];
% 
% k=1;
% for i=1:length(foldernames)
%     disp([i/length(foldernames)]);
%     folderpath = fullfile(root_path, foldernames{i},'GangulyServer',foldernames{i},'Robot3DArrow');
%     D=dir(folderpath);
%     if length(D)==0
%         folderpath = fullfile(root_path, foldernames{i},'GangulyServer','Robot3DArrow');
%         D=dir(folderpath);
%     end
%     task_files_temp=[];
%     for j=3:length(D)
%         filepath=fullfile(folderpath,D(j).name,'BCI_Fixed');
%         if exist(filepath)
%             task_files_temp = [task_files_temp;findfiles('mat',filepath)'];
%         end
%     end
%     if ~isempty(task_files_temp)
%         task_files = [task_files;task_files_temp];k=k+1;
%     end
%     %     imag_files{i} = imag_files_temp;
%     %     online_files{i} = online_files_temp;
% end


clearvars -except task_files
close all

files_not_loaded=[];
% good_dis=[];bad_dis=[];
for file = 501:2000

    if mod(file, 50) == 0
        disp([file/length(task_files)]);
    end

    try
        load(task_files{file})
        file_loaded = true;
    catch
        file_loaded=false;
        disp(['Could not load ' task_files{file}]);
        files_not_loaded=[files_not_loaded;task_files(file)];
    end
    if file_loaded

        result=TrialData.ClickerState==TrialData.TargetID;
        decoding_acc=length(find(result==1))/length(result);

        if decoding_acc>0.3 &&decoding_acc<0.8
            continue
        end

        raw = cell2mat(TrialData.BroadbandData');
        chMap = TrialData.Params.ChMap;
        Fs = 1000; % 1000 Hz

        bpFilt = designfilt('bandpassiir','FilterOrder',4, ...
            'HalfPowerFrequency1',3,'HalfPowerFrequency2',40, ...
            'SampleRate',1e3);

        peakFreqs = zeros(1, 128);
        for i=1:128
            data = raw(:,i);

            [wt,f] = cwt(data', 'amor', Fs); % 'amor' - Morlet

            power = mean((abs(wt).^2),2);


            logf = log10(f);
            logPower = log10(power);

            % fit
            [bhat p wh se ci t_stat]=robust_fit(logf,logPower,1);


            % figure
            % plot(logf,logPower)
            % hold on; plot(logf,bhat(1)+bhat(2).*logf,'r')
            %
            % figure;
            % plot(logf,logPower - bhat(1)+bhat(2).*logf)

            % Detrended log power
            detrendedLogPower = logPower - (bhat(1) + bhat(2).*logf);

            % Find peaks
            [pks, locs] = findpeaks(detrendedLogPower, 'MinPeakProminence', 0.1);

            % Check if peaks were found
            if ~isempty(pks)
                % Find the frequency of the largest peak
                [~, maxIdx] = max(pks);
                peakFreqs(i) = round(f(locs(maxIdx)));
            else
                peakFreqs(i) = NaN;
            end
        end

        bpFilt = designfilt('bandpassiir','FilterOrder',4, ...
            'HalfPowerFrequency1',mode(peakFreqs)-1,'HalfPowerFrequency2',mode(peakFreqs)+1, ...
            'SampleRate',1e3);

        filteredSignal = filtfilt(bpFilt,raw);

        filteredSignal = filteredSignal';

        phaseMatrix = zeros(size(filteredSignal));
        for i = 1:size(filteredSignal, 1)
            tmp = filteredSignal(i, :);
            phaseMatrix(i, :) = angle(hilbert(tmp));
        end


        PLV = zeros(1, size(phaseMatrix,2));
        for t = 1:size(phaseMatrix,2)
            phase_differences = exp(1i * (phaseMatrix(:, t) - mean(phaseMatrix(:, t))));
            PLV(t) = abs(mean(phase_differences));
        end
        tmp = find(TrialData.TaskState==3);
        segments=[];
        for i=1:length(tmp)
            tmp_start = size(cell2mat(TrialData.BroadbandData(1:tmp(i)-1)'),1);
            tmp_end = size(cell2mat(TrialData.BroadbandData(1:tmp(i))'),1);
            tmp_segments = [tmp_start, tmp_end, TrialData.ClickerState(i)==TrialData.TargetID];
            segments = [segments;tmp_segments];
        end

        figure;
        time = 1:size(phaseMatrix,2);
        plot(time, PLV);
        xlabel('Time');
        ylabel('Phase Locking Value across Channel');
        title('PLV for coherence over Time');


        hold on;

        for i = 1:size(segments, 1)
            startPt = segments(i, 1);
            endPt = segments(i, 2);
            colorCode = segments(i, 3);

            if colorCode == 0
                color = 'red';
            else
                color = 'green';
            end

            patch([startPt, endPt, endPt, startPt], [0, 0, max(PLV), max(PLV)], color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end

        hold off;

        cd('F:\BCICodeGroup\NikTravelingWave\fig');

        [peaks_max, locations_max] = findpeaks(PLV);
        [peaks_min, locations_min] = findpeaks(-PLV);
        locations = [locations_max, locations_min];
        kinax = TrialData.TaskState; state1 = find(kinax==1);
        state2 = find(kinax==2); state3 = find(kinax==3);
        state4 = find(kinax==4);

        len_state1 = size(cell2mat(TrialData.BroadbandData(state1)'),1);
        len_state2 = size(cell2mat(TrialData.BroadbandData(state2)'),1)+len_state1;
        len_state3 = size(cell2mat(TrialData.BroadbandData(state3)'),1)+len_state2;

        distances = abs(locations - len_state2);
        [minDistance, index] = min(distances);

        if decoding_acc<0.3
            filename = sprintf('bad_trial_%d.png', file);
            bad_dis = [bad_dis minDistance];
        else
            filename = sprintf('good_trial_%d.png', file);
            good_dis = [good_dis minDistance];
        end


        saveas(gcf, filename);
        cd('F:\BCICodeGroup\NikTravelingWave');
        close(gcf);
    end

end




