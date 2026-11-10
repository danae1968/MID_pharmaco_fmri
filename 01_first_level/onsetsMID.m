function [onsets, durations, trialMat100] = onsetsMID
% ONSETSMID  Extract MID task event onsets from Presentation log files.
%
% For every subject and drug session, reads the Presentation .log file of
% the monetary incentive delay (MID) task, finds the cue, target and
% feedback events, and saves an SPM "multiple conditions" file
% (names/onsets/durations) next to the log file:
%
%   bids/sub-XXX/ses-drugN/func/sub-XXX_ses-drugN_task-MID_onsets.mat
%
% Conditions saved (all modelled as events, duration 0):
%   LowCue, HighCue   - reward cue (low / high incentive)
%   Target            - target onset
%   LowHit, HighHit   - feedback after a hit (low / high incentive)
%   LowMiss, HighMiss - feedback after a miss (low / high incentive)
%   Resp              - button presses (sanity-check regressor)
%
% Onsets are in seconds relative to the first trial after the dummy scans.
% This file is the input for the 'multipleConditions' step of
% firstLevel_MID_fir.m.

% Directories
% -------------------------------------------------------------------------
if ispc
    projectDir = 'P:\3017048.01\';
elseif isunix
    projectDir = '/project/3017048.01/';
else
    error('Unknown OS')
end

trialMat100 = [];
bidsDir = fullfile(projectDir,'bids');

% Subjects and sessions
% -------------------------------------------------------------------------
subNo = 1:100;
ses   = repmat([{{1 2 3}}],length(subNo),1);

saveF    = 1;     % write the onsets .mat file
prescanN = 30;    % Presentation events logged before the first real trial

for ii = 1:numel(subNo)
    i = subNo(ii);
    for jj = 1:numel(ses{ii})
        j = ses{ii}{jj};

        subDir    = fullfile(bidsDir,sprintf('sub-%.3d',i),sprintf('ses-drug%.d',j),'func');
        filename  = dir(fullfile(subDir,sprintf('*sub-%.3d_ses-drug%d_task-MID_*.log',i,j)));
        outputDir = subDir;

        fid = fopen(fullfile(subDir,filename.name),'r');

        if fid ~= - 1

            skip = 3;

            % Read the tab-delimited Presentation log: header row, then data
            hdr  = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s', 1, 'delimiter', '\t');
            data = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s','delimiter', '\t');

            % Event codes as written by the Presentation script
            low_cue   = 'low_mon_rw_cue';
            high_cue  = 'high_mon_rw_cue';
            target    = '_tar_';
            hit       = 'hit';
            miss      = 'miss';
            response  = 'Response';
            low_hit   = 'low_rw_hit';
            high_hit  = 'high_rw_hit';
            high_miss = 'high_rw_mis';
            low_miss  = 'low_rw_mis';

            % The experimenter presses '5' at the end of the practice block;
            % the first real trial follows after the dummy-scan events.
            endprac_idx = find(strcmp(data{4}, '5'));
            startRaw = endprac_idx+prescanN+1;

            Otime    = str2double(data{5}(startRaw:end));    % event times
            tStart   = Otime(1);
            Pic      = data{4}(startRaw:end);                % event code
            Duration = data{8}(startRaw:end);
            Event    = data{3}(startRaw:end);                % event type
            Code     = zeros(length(Pic),1);                 % numeric event code (below)
            Resp     = zeros(length(Event),1);               % 1 = button press
            eventNumber = str2double(erase([data{2}(startRaw:end)],'A'));   % trial number

            for w = 1:length(Event)
                if strcmp(Event{w},response)
                    Resp(w) = 1;
                end
            end

            % Numeric code per event:
            %   1 low cue, 2 high cue, 3 target, 5 low hit, 6 high hit,
            %   7 hit, 8 button press, 9 high miss, 10 low miss
            for ind = 1:length(Pic)
                if ischar(Pic{ind})
                    if strcmp(Pic{ind},low_cue)
                        Code(ind) = 1;
                    elseif strcmp(Pic{ind},high_cue)
                        Code(ind) = 2;
                    elseif regexpi(Pic{ind},target)
                        Code(ind) = 3;
                    elseif strcmp(Pic{ind},low_hit)
                        Code(ind) = 5;
                    elseif strcmp(Pic{ind},high_hit)
                        Code(ind) = 6;
                    elseif strcmp(Pic{ind},hit)
                        Code(ind) = 7;
                    elseif strcmp(Pic{ind},'1')
                        Code(ind) = 8;
                    elseif strcmp(Pic{ind},high_miss)
                        Code(ind) = 9;
                    elseif strcmp(Pic{ind},low_miss)
                        Code(ind) = 10;
                    end
                end
            end

            % Presentation logs time in units of 0.1 ms -> divide by 10000
            % for seconds.
            OnsetsLH = (Otime(Code==5)-tStart)/10000;
            OnsetsHH = (Otime(Code==6)-tStart)/10000;
            DurHH    = cell2mat(Duration(Code==6))/10000;
            OnsetsHM = (Otime(Code==9)-tStart)/10000;
            OnsetsLM = (Otime(Code==10)-tStart)/10000;
            DurHM    = cell2mat(Duration(Code==9))/10000;
            OnsetsLC = (Otime(Code==1)-tStart)/10000;
            DurLC    = cell2mat(Duration(Code==1))/10000;
            OnsetsHC = (Otime(Code==2)-tStart)/10000;
            DurHC    = cell2mat(Duration(Code==2))/10000;
            OnsetsC  = (Otime(Code==2 | Code==1)-tStart)/10000;
            OnsetsT  = (Otime(Code==3)-tStart)/10000;
            DurT     = cell2mat(Duration(Code==3))/10000;
            OnsetsH  = (Otime(Code==7)-tStart)/10000;
            OnsetsR  = (Otime(Resp==1)-tStart)/10000;
            OnsetsR  = (Otime(Code==8)-tStart)/10000;

            trialNum    = length(OnsetsT);
            trialOnsets = sort([OnsetsHC;OnsetsLC]);
            hits        = sum(length(OnsetsH));
            trialEvNum  = eventNumber(Code==3);
            trialType   = Code(Code==1 | Code==2, :);   % 1 = low reward, 2 = high reward
        end

        % Classify each trial's response
        %   -1 premature (before target), 1 hit, 2 too late, 0 no response
        for trial = 1:trialNum

            if fid ~= -1
                % Timings of all events in this trial
                trialTimings = (Otime(eventNumber==trialEvNum(trial))-tStart)/10000;
                trialDur = trialTimings(end)-trialTimings(1);

                if sum(Resp(eventNumber==trialEvNum(trial)))>0
                    trialRT(trial)    = trialTimings(find(Resp(eventNumber==trialEvNum(trial)),1))-OnsetsT(trial);
                    trialRTCue(trial) = trialTimings(find(Resp(eventNumber==trialEvNum(trial)),1))-OnsetsC(trial);
                    if trialRT(trial)<0
                        trialResp(trial) = -1;
                    elseif sum(Code(eventNumber==trialEvNum(trial))==7)>0
                        trialResp(trial) = 1;
                    else
                        trialResp(trial) = 2;
                    end

                elseif sum(Resp(eventNumber==trialEvNum(trial)))==0
                    trialResp(trial)  = 0;
                    trialRT(trial)    = NaN;
                    trialRTCue(trial) = NaN;
                end

            elseif fid==-1
                warning(sprintf('filename does not exist for subject %d session %d',i,j))
                trialResp(trial)  = NaN;
                trialRT(trial)    = NaN;
                trialRTCue(trial) = NaN;
                trialType(trial)  = NaN;
            end
        end

        if fid>0
            fclose(fid);
        end

        % Consistency check: classified hits vs hit feedback events
        if sum(trialResp==1)~=hits
            warning(sprintf('The number of correct trials does not match the hits for participant %d session %d',i,j))
        end

        % Cue onsets split by trial outcome (not used in the current model)
        OnsetsHitC  = OnsetsC(trialResp==1);
        OnsetsHHC   = intersect(OnsetsHC,OnsetsHitC);
        OnsetsLHC   = intersect(OnsetsLC,OnsetsHitC);
        OnsetsMissC = OnsetsC(trialResp~=1);
        OnsetsHMC   = intersect(OnsetsHC,OnsetsMissC);
        OnsetsLMC   = intersect(OnsetsLC,OnsetsMissC);

        % Save SPM multiple-conditions file
        names     = {'LowCue' 'HighCue' 'Target' 'LowHit' 'HighHit' 'LowMiss' 'HighMiss' 'Resp'};
        onsets    = {OnsetsLC OnsetsHC OnsetsT OnsetsLH OnsetsHH OnsetsLM OnsetsHM OnsetsR};
        durations = {0 0 0 0 0 0 0 0};
        outputfileRewAnt = fullfile(outputDir,sprintf('sub-%.3d_ses-drug%d_task-MID_onsets.mat',i,j));

        if saveF
            save(outputfileRewAnt,'names','onsets','durations')
        end
    end
end

end
