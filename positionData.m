classdef positionData < ephysData
    %% Class to perform and consolidate analysis relating to bat spatial position data
    % Initializing this class will collect all available tracking data and
    % organize it according to experiment date (expDate) and bat ID number
    % (batNum). Then: 1) raw LED positions are read in 2) corresponding
    % video timestamps are read in 3) positions are organized by bat, and
    % 4) cleaned up (centered, gaps filled, etc).
    %
    % The main data structures in a pData object are 'batPos' and 'posTS.'
    % Both are MATLAB "container.Maps" objects (i.e. dictionaries) that can
    % be accessed by unique keys. posTS is indexed by expDate and returns
    % the timestamps of the corresponding positions for all bats for that
    % day. batPos is nested, with the first level of nesting indexed by
    % expDate, and the next level by batNum.
    %%
    properties
        gap_fill_window_s = 2 % how long (in s) to fill in gaps in LED tracking
        video_fs = 20 % sampling rate of video cameras
        all_bat_nums = [11636,11682,13688,14612,14654,14798,60064,71216,71335] % all bat numbers we'd like to keep track of
        sessionType = 'social' % only 'social' sessions work for now
        callOffset = 2 % time offset (in s) +/- around calls over which to average bat position
        pixel2cm = 0.25 % factor to convert pixels to cm
        
        batPos % XY positions of each bat in a nested map, indexed by expDate and then batNum
        posTS % timestamps of positions in a map indexed by expDate
        foodTime
        tracking_data_dir
        video_data_dir
        call_data_dir
        
    end
    
    methods
        % main utility functions
        function pd = positionData(varargin)
            %% Initializes positionData
            % Inputs:
            % used_exp_dates: list of datetimes to limit analysis to
            % used_bat_nums: list of bat IDs to limit analysis to
            
            pnames = {'used_exp_dates','used_bat_nums'};
            dflts  = {[],[]};
            [used_exp_dates,used_bat_nums] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            pd = pd@ephysData('adult_social'); % base pd on ephysData to get basic experimental metadata
            pd.tracking_data_dir = fullfile(pd.baseDirs{1},'tracking_data');
            pd.video_data_dir = fullfile(pd.baseDirs{1},'video_data');
            
            mov_window_samples = pd.gap_fill_window_s*pd.video_fs; % get number of frames over which to perform moving window cleaning of positions
            
            % get a list of all LEDtracking results and the corresponding
            % expDates
            led_track_fnames = dir(fullfile(pd.tracking_data_dir,['LEDtracking_pred_' pd.sessionType '*.mat']));
            fnameSplit = arrayfun(@(x) strsplit(x.name,'_'),led_track_fnames,'un',0);
            exp_date_strs = cellfun(@(x) x{end}(1:end-4),fnameSplit,'un',0);
            expDates = pd.expstr2datetime(exp_date_strs);
            
            if ~isempty(used_exp_dates) % if supplied, limit expDates to those supplied
                [~,used_date_idx] = ismember(expDates,used_exp_dates);
                used_date_idx = logical(used_date_idx);
                expDates = expDates(used_date_idx);
                exp_date_strs = exp_date_strs(used_date_idx);
                led_track_fnames = led_track_fnames(used_date_idx);
            end
            
            if ~isempty(used_bat_nums) % if supplied, limit batNums to those supplied
                pd.all_bat_nums = used_bat_nums;
            end
            
            nExp = length(exp_date_strs);
            used_exp_idx = true(1,nExp);
            [bat_pos_cell, pos_ts_cell] = deal(cell(1,nExp));
            foodTimes = nan(1,nExp);
            
            for exp_k = 1:nExp
                led_track_fname = fullfile(led_track_fnames(exp_k).folder,led_track_fnames(exp_k).name);
                exp_date_str = exp_date_strs{exp_k};
                % get 'frame_ts_info' file which contains timestamps for
                % each video frame, if this file doesn't exist, skip this
                % date
                frame_ts_info_fname = fullfile(pd.video_data_dir,[exp_date_str '_color_frame_timestamps_info_social.mat']);
                
                if ~isfile(frame_ts_info_fname)
                    used_exp_idx(exp_k) = false;
                    continue
                end
                
                % get food session time
                foodTimes(exp_k) = get_food_time(pd,expDates(exp_k));
                
                % load frame_ts_info and LEDtracks
                s = load(frame_ts_info_fname);
                frame_ts_info = s.frame_ts_info;
                LEDTracks = load(led_track_fname,'file_frame_number','fileIdx',...
                    'color_pred_model','centroidLocs','predPosterior','predColors');
                
                % get index of frames for which we have timestamps and LED
                % tracking results(should be ~99.9% of frames)
                [idx_tracks,idx_ts] = align_tracks_with_frame_ts(LEDTracks,frame_ts_info);
                
                % only use timestamps that overlap with the frames
                % processed in LED tracking
                timestampsNlg = frame_ts_info.timestamps_nlg(idx_ts);
                
                % make sure the provided sampling rate is correct
                assert(round(1e3/median(diff(timestampsNlg))) == pd.video_fs)
                
                % get a matrix of XY positions for each frame and bat and
                % deal with duplicate colors within frames
                predCentroids = get_pred_centroids(LEDTracks);
                % only use LED tracks that overlap with the timestamps
                predCentroids = predCentroids(idx_tracks,:,:);
                % filter, center, and fill gaps in LED tracks
                predCentroids = clean_pred_centroids(predCentroids,'GapMethod','movmedian','movWindow',mov_window_samples);
                
                % impose the same order on the position matrix for each day
                colorStrs = LEDTracks.color_pred_model.ClassificationSVM.ClassNames;
                predCentroids = reorder_bat_pos(pd,predCentroids,colorStrs,expDates(exp_k));
                % convert pixel values to cm
                predCentroids = pd.pixel2cm*predCentroids;
                % convert to cell before making into a Map
                predCentroids = squeeze(num2cell(predCentroids,[1 2]));
                % each expDate's batPos Map is indexed by batNum (a double,
                % not a string)
                bat_pos_cell{exp_k} = containers.Map(pd.all_bat_nums,predCentroids);
                pos_ts_cell{exp_k} = timestampsNlg';
                
            end
            
            % both batPos and posTS are indexed by expDate strings of the
            % form 'mmddyyyy'
            pd.batPos = containers.Map(exp_date_strs(used_exp_idx),bat_pos_cell(used_exp_idx));
            pd.posTS = containers.Map(exp_date_strs(used_exp_idx),pos_ts_cell(used_exp_idx));
            pd.foodTime = containers.Map(exp_date_strs(used_exp_idx),foodTimes(used_exp_idx));
            
        end
        function pos = get_pos(pd,expDate,batNum,varargin)
            % utility function to get bat position(s) for a given expDate
            % and, optionally, for a given batNum. expDate can be either a
            % datetime or a date string. Without a batNum, this returns
            % the map of bat positions for the expDate, with a batNum this
            % returns the array of positions for that bat on that date
            
            pnames = {'sessionSelection'};
            dflts  = {'social'};
            [sessionSelection] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            if isdatetime(expDate)
                expDate = datestr(expDate,'mmddyyyy');
            end
            pos = pd.batPos(expDate);
            if nargin > 2
                if ~isnumeric(batNum)
                    batNum = str2double(batNum);
                end
                pos = pos(batNum);
                fTime = pd.foodTime(expDate);
                t = pd.posTS(expDate);
                switch sessionSelection
                    case 'social'
                        if isnan(fTime)
                            return
                        end
                        idx = t < fTime;
                    case 'food'
                        if isnan(fTime)
                            pos = [];
                            return
                        end
                        idx = t > fTime;
                    case 'all'
                        idx = true(1,length(t));
                end
                pos = pos(idx,:);
            end
        end
        function t = get_time(pd,expDate)
            % utlity function get timestamps for a given expDate. expDate
            % can be either a datetime or a date string.
            if isdatetime(expDate)
                expDate = datestr(expDate,'mmddyyyy');
            end
            t = pd.posTS(expDate);
        end
        function dist = get_dist(pd,batNums,expDate,varargin)
            % gets the frame by frame euclidean distance between two bats
            % on a given expDate
            % Inputs:
            % batNums: 1 x 2 array of batNums
            % expDate: datetime or date string
            
            pnames = {'sessionSelection'};
            dflts  = {'social'};
            [sessionSelection] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            dist = vecnorm((get_pos(pd,expDate,batNums(1),'sessionSelection',sessionSelection)...
                - get_pos(pd,expDate,batNums(2),'sessionSelection',sessionSelection)),2,2);
        end
        function pairDist = get_pairwise_dist(pd,expDate,varargin)
            % gets pairwise distance between all pairs of bats.
            % Inputs:
            % used_bat_nums: optional list to limit which bats to look at
            % Outputs:
            % pairDist: a new map indexed by strings of the form
            % '[batNum1]-[batNum2]' whose values are the frame by frame
            % distance between that pair.
            
            pnames = {'used_bat_nums','sessionSelection'};
            dflts  = {[],'social'};
            [used_bat_nums,sessionSelection] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            % first enumerate all possible bat pairs
            batPairs = get_bat_pairs(pd,'expDate',expDate,'used_bat_nums',used_bat_nums);
            
            % get the distance between each of those pairs
            nPairs = size(batPairs,1);
            pairDist = cell(1,nPairs);
            for bat_pair_k = 1:nPairs
                pairDist{bat_pair_k} = get_dist(pd,batPairs(bat_pair_k,:),expDate,'sessionSelection',sessionSelection);
            end
            
            % convert the array of bat pairs into a list of strings of the
            % form '[batNum1]-[batNum2]' to use as keys into the Map of
            % pairwise distances
            bat_pair_keys = pd.get_pair_keys(batPairs);
            pairDist = containers.Map(bat_pair_keys,pairDist);
            
        end
        function batPairs = get_bat_pairs(pd,varargin)
            % utility function to enumerate all bat pairs
            % Inputs:
            % expDate: optional datetime or date string to only list pairs
            % of bats that were present on this day
            % used_bat_nums: optional list to limit which bats to look at
            % Outputs:
            % batPairs: (nBat choose 2) x 2 array of bat numbers listing
            % all possible bat pairs
            pnames = {'expDate','used_bat_nums'};
            dflts  = {[],[]};
            [expDate,bat_num_list] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            used_bat_nums = pd.all_bat_nums;
            
            % if expDate is provided get a list of batNums used on this
            % date and limit bats to just those
            if ~isempty(expDate)
                bat_color_table = get_bat_color_table(pd,expDate);
                used_bat_nums = intersect(used_bat_nums,bat_color_table.batNum);
            end
            
            if ~isempty(bat_num_list)
                used_bat_nums = intersect(used_bat_nums,bat_num_list);
            end
            
            batPairs = nchoosek(used_bat_nums,2);
            
        end
        function bat_pair_keys = get_pair_keys(~,batPairs)
            % utility function to convert a nPair x 2 array of batNums to a
            % list of strings of the form '[batNum1]-[batNum2]'
            bat_pair_keys = cellfun(@(batPair) num2str(sort(batPair),'%d-%d'),num2cell(batPairs,2),'UniformOutput',false);
        end
        function batPairs = get_key_pairs(~,bat_pair_keys)
            % utility function to convert a list of bat pair strings back
            % to a nPair x 2 array of batNums
            batPairs = cellfun(@(pairKey) str2double(strsplit(pairKey,'-')),bat_pair_keys,'un',0);
            batPairs = vertcat(batPairs{:});
        end
        function exp_date_strs = datetime2expstr(~,expDates)
            % utility function to convert datetimes to date strings.
            % returns as cell array of strings, unless only one date is
            % provided.
            exp_date_strs = cellfun(@(expDate) datestr(expDate,'mmddyyyy'),num2cell(expDates),'un',0);
            if length(exp_date_strs) == 1
                exp_date_strs = exp_date_strs{1};
            end
        end
        function expDates = expstr2datetime(~,exp_date_strs)
            % utility function to convert date strings into datetimes
            if ~iscell(exp_date_strs)
                exp_date_strs = {exp_date_strs};
            end
            expDates = cellfun(@(dateStr) datetime(dateStr,'InputFormat','MMddyyyy'),exp_date_strs);
        end
        
        % functions to get call related position/distance data
        function callPos = get_call_pos(pd,cData,varargin)
            % gets the position of bats averaged around the time that a
            % call occurred across expDates.
            % Inputs:
            % cData: callData object corresponding to this experiment
            % inter_call_int: minimum interval (in ms) between call
            % occurrences, if provided, all calls separated by less than
            % that amount are excluded
            % expDates: if provided, only look at these expDates. defaults
            % to using all expDates.
            % Outputs:
            % callPos: nested map, first indexed by expDate, then by
            % batNum. Values are structs with fields 'pos' and 'caller'
            % which contain the given bat's position and which bat made the
            % call, respectively
            
            pnames = {'inter_call_int','expDates','call_offset_range','avgPos'};
            dflts  = {-Inf,[],[],true};
            [inter_call_int,expDates,call_offset_range,avgPos] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            % set up +/- offset in ms
            if isempty(call_offset_range)
                call_offset_range = pd.callOffset.*[-1 1];
            end
            
            call_sample_offset = round(call_offset_range*pd.video_fs);
            call_sample_offset = call_sample_offset(1):call_sample_offset(2)-1;
            
            % get expDates in datetime and corresponding date strings
            if isempty(expDates)
                exp_date_strs = pd.batPos.keys;
                expDates = pd.expstr2datetime(exp_date_strs);
            else
                exp_date_strs = pd.datetime2expstr(expDates);
            end
            
            % iterate over all expDates
            nExp = length(expDates);
            callPos = cell(1,nExp);
            for exp_k = 1:nExp
                pos = pd.get_pos(expDates(exp_k)); % get position of all bats on this expDate
                callTimes = cData('expDay',expDates(exp_k)).callPos; % returns array of call start and stop times (in s) that occurred on this expDate
                callTimes = 1e3*callTimes(:,1); % use only call start times and convert to ms
                calling_bat_nums = cData('expDay',expDates(exp_k)).batNum'; % gets the batID that produced these calls
                callIDs = cData('expDay',expDates(exp_k)).callID';
                
                % if provided, limit to calls separated by minimum inter
                % call interval (default is -Inf)
                callIdx = [Inf; diff(callTimes)] > inter_call_int;
                callTimes = callTimes(callIdx);
                calling_bat_nums = calling_bat_nums(callIdx);
                callIDs = callIDs(callIdx);
                
                % get the frame timestamps for this expDate
                t = pd.get_time(expDates(exp_k));
                
                % here we'll iterate over all calls and then over bats
                nCall = length(callTimes);
                nBat = length(pd.all_bat_nums);
                call_pos_idx = cell(1,nCall);
                current_call_pos = cell(1,nBat);
                
                % for each call get the corresponding index into the
                % position data over which to average
                for call_k = 1:nCall
                    [~,current_call_sample] = min(abs(t - callTimes(call_k)));
                    call_pos_idx{call_k} = current_call_sample + call_sample_offset;
                    if any(call_pos_idx{call_k} < 0 | call_pos_idx{call_k} > length(t))
                       call_pos_idx{call_k} = []; 
                    end
                end
                
                % for each bat, get its average position for all calls and
                % store as a struct with fields 'pos' and 'caller' which
                % contain the given bat's position and which bat made the
                % call (we might want to save also the unique callID here)
                for bat_k = 1:nBat
                    current_bat_pos = pos(pd.all_bat_nums(bat_k));
                    if avgPos
                        current_bat_pos = cellfun(@(idx) nanmean(current_bat_pos(idx,:)),call_pos_idx,'un',0);
                    else
                        current_bat_pos = cellfun(@(idx) current_bat_pos(idx,:),call_pos_idx,'un',0);
                    end
                    [current_bat_pos{cellfun(@isempty,current_bat_pos)}] = deal(nan(length(call_sample_offset),2));
                    current_call_pos{bat_k} = struct('pos',current_bat_pos,'caller',calling_bat_nums,'callID',num2cell(callIDs));
                end
                % create map of structs indexed by batNums
                callPos{exp_k} = containers.Map(pd.all_bat_nums,current_call_pos);
            end
            % create a map of maps indexed by expDate
            callPos = containers.Map(exp_date_strs,callPos);
        end
        function callDist = get_call_dist(pd,cData,varargin)
            % gets the distance between pairs of bats around the time that
            % a call occurred across expDates.
            % Inputs:
            % cData: callData object corresponding to this experiment
            % inter_call_int: minimum interval (in ms) between call
            % occurrences, if provided, all calls separated by less than
            % that amount are excluded
            % expDates: if provided, only look at these expDates. defaults
            % to using all expDates.
            % Outputs:
            % callDist: nested map, first indexed by expDate, then by
            % bat pair string. Values are structs with fields 'dist' and
            % 'caller' which contain the given bat pair's distance and
            % which bat made the call, respectively
             pnames = {'inter_call_int','expDates','call_offset_range','avgPos'};
            dflts  = {-Inf,[],[],true};
            [inter_call_int,expDates,call_offset_range,avgPos] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            % get expDates in datetime and corresponding date strings
            if isempty(expDates)
                exp_date_strs = pd.batPos.keys;
                expDates = pd.expstr2datetime(exp_date_strs);
            else
                exp_date_strs = pd.datetime2expstr(expDates);
            end
            
            % get bat positions around calls across expDates
            call_bat_pos = get_call_pos(pd,cData,'expDate',expDates,'inter_call_int',inter_call_int,'call_offset_range',call_offset_range,'avgPos',avgPos);
            
            nExp = length(expDates);
            callDist = cell(1,nExp);
            for exp_k = 1:nExp
                % enumerate list of bat pairs present on this day
                batPairs = get_bat_pairs(pd,'expDate',expDates(exp_k));
                nPairs = size(batPairs,1);
                % get the map of bat position indexed by batNum for this
                % expDate
                current_bat_pos = call_bat_pos(exp_date_strs{exp_k});
                bat_call_dist = cell(1,nPairs);
                % iterate over all bat pairs
                for bat_pair_k = 1:nPairs
                    % get the struct of position and calling batNum of both
                    % bats in this bat pair
                    current_call_bat_pos = cellfun(@(bat) current_bat_pos(bat),num2cell(batPairs(bat_pair_k,:)),'un',0);
                    
                    % get the calling bat ID for each call and assert that
                    % the list of calling bats is the same across bats
                    callIDs = cellfun(@(bat) [bat.callID],current_call_bat_pos,'un',0);
                    assert(isequal(callIDs{:}))
                    
                    calling_bat_nums = cellfun(@(bat) {bat.caller},current_call_bat_pos,'un',0);
                    % get the array of positions for both bat in this pair
                    current_call_bat_pos = cellfun(@(bat) cat(3,bat.pos),current_call_bat_pos,'un',0);
                    % calculate the distance between this pair of bats
                    current_call_dist = squeeze(vecnorm(current_call_bat_pos{1} - current_call_bat_pos{2},2,2));
                    % save as a struct with fields 'dist' and 'caller'
                    bat_call_dist{bat_pair_k} = struct('dist',num2cell(current_call_dist,1),'caller',calling_bat_nums{1},'callID',num2cell(callIDs{1}));
                end
                % create map of structs indexed by bat pair strings
                bat_pair_keys = pd.get_pair_keys(batPairs);
                callDist{exp_k} = containers.Map(bat_pair_keys,bat_call_dist);
            end
            % create a map of maps indexed by expDate
            callDist = containers.Map(exp_date_strs,callDist);
        end
        function callMap = collect_by_calls(~,call_data_map)
            %% reorganizes a call-distance/position map from indexing by expDates to
            % indexing by callIDs
            % Inputs:
            % call_data_map: output of pd.get_call_dist or pd.get_call_pos
            % Outputs:
            % callMap: nested map indexed by callID and batPair which
            % contains pairwise bat distances/positions
            
            % first determine the variable name of the data we're working
            % with ('pos' or 'dist')
            firstMap = call_data_map.values;
            firstMap = firstMap{1};
            firstMap = firstMap.values;
            firstMap = firstMap{1};
            varName = setdiff(fieldnames(firstMap),{'caller','callID'});
            varName = varName{1};
            
            callMap = containers.Map('KeyType','double','ValueType','any');
            % iterate over expDates
            for expKey = call_data_map.keys
                expData = call_data_map(expKey{1});
                % iterate over bat pairs in this expDate
                for batKey = expData.keys
                    batData = expData(batKey{1});
                    callIDs = [batData.callID];
                    call_k = 1;
                    % iterate separately over callIDs for each batPair
                    for call_id = callIDs
                        if isKey(callMap,call_id) % if we've already added this callID for a different batPair
                            bat_call_map = callMap(call_id); % get this call's distance map
                            bat_call_map(batKey{1}) = batData(call_k).(varName); % get this bat pair's distance for this call
                            callMap(call_id) = bat_call_map; % replace this call's distance map with updated map
                        else % if we haven't already, initialize a map for this call to be indexed by batPair
                            callMap(call_id) = containers.Map(batKey{1},batData(call_k).(varName));
                        end
                        call_k = call_k + 1;
                    end
                end
            end
        end
        
        % functions to analyze distance alone
        function [dist_matrix_corr, mantel_pval] = compare_dist_matrices(pd,varargin)
            %% gets correlation values between distance matrices across days and corresponding Mantel Test p-values 
            % Inputs:
            % nRep: number of permutations to perform for p-values
            % exclDates: dates to exclude (NOTE: this analysis requires
            % that all dates have the same bats)
            % Outputs:
            % dist_matrix_corr: correlation values between each expDates
            % distance matrix and the average distance matrix across all
            % other dates
            % mantel_pval: corresponding p-values for those correlation
            pnames = {'nRep','exclDates'};
            dflts  = {1e4,'08032020'};
            [nRep,exclDates] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            exp_date_strs = pd.batPos.keys;
            exp_date_strs = setdiff(exp_date_strs,exclDates);
            expDist = cellfun(@(expDay) posData.get_pairwise_dist(expDay),exp_date_strs,'un',0);
            
            assert(all(cellfun(@(x) all(strcmp(x.keys,expDist{1}.keys)),expDist)))
            
            expDist = cellfun(@(dist) squareform(cellfun(@nanmean,dist.values)),expDist,'un',0);
            [dist_matrix_corr, mantel_pval]= cellfun(@(m,k) bramila_mantel(m,nanmedian(cat(3,expDist{setdiff(1:length(expDist),k)}),3),nRep,'pearson'),expDist,num2cell(1:length(expDist)));
        end
        
        % functions to analyze call distance and inter-brain correlation
        function [predTable,call_dist_and_corr] = get_call_dist_and_corr(pd,cData,bat_pair_corr_info,call_pair_type,varargin)
            %% gets pairwise distance and interbrain correlation between bats during calls
            % NOTE: This function is meant to analyze the relationship
            % between inter-brain correlation around calls as a function of
            % pairwise distance. For caller-to-listener ('c2l') correlation,
            % this can be modeled as a straightforward correlation between
            % distance and correlation between all c2l pairs:
            % R(CL) ~ D(CL)
            % For listener-to-listener ('l2l') correlation, the distance
            % between listeners also needs to be taken into account:
            % R(L1L2) ~ D(L1L2) + D(L1C) + D(L2C)
            % Inputs:
            % cData: callData object
            % bat_pair_corr_info: output of 'calculate_all_cross_brain_lfp_corr'
            % which contains the instantaneous inter-brain correlation 
            % between bats around calls.
            % call_pair_type: Specifies  which type of inter-brain 
            % correlation to look at. Either 'l2l' (listener-to-listener) 
            % or 'c2l' (caller-to-listener). 
            % inter_call_int: sets which calls to look at by
            % inter-call interval
            % tLim: time limits relative to call onset over which to
            % average inter-brain correlation
            % Ouputs: 
            % predTable: a table that can input into 'fitlm' to test the
            % relationship between distance and interbrain correlation.
            % call_dist_and_corr: a structure containing maps indexed by
            % callID listing the correlation and distance between bat pairs
            
            pnames = {'inter_call_int','tLim'};
            dflts  = {3e3,[-0.3 0.3]};
            [inter_call_int,tLim] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            % get a map of call distance, indexed by callID
            callDist = pd.get_call_dist(cData,'inter_call_int',inter_call_int);
            call_dist_by_calls = pd.collect_by_calls(callDist);
            
            % convert the bat_pair_corr_info struct into a map indexed by
            % callID. This also returns a map of all the callIDs included 
            % in each call bout that weren't used explicitly because their
            % inter_call_interval was too short
            [bat_pair_corr_map, included_call_map] = bat_pair_corr_info_2_map(bat_pair_corr_info,pd);
            
            % get a list of callIDs to iterate over
            all_call_ids = call_dist_by_calls.keys;
            all_call_ids = [all_call_ids{:}];
            
            % get the index into the correlation around calls by time
            [~,t_idx] = inRange(bat_pair_corr_info.time,tLim);
            k = 1;
            
            % initialize maps
            l2lCorr = containers.Map('KeyType','double','ValueType','any');
            c2lDist = containers.Map('KeyType','double','ValueType','any');
            l2lDist = containers.Map('KeyType','double','ValueType','any');
            c2lCorr = containers.Map('KeyType','double','ValueType','any');
            
            % iterate over call IDs
            for callID = all_call_ids
                if isKey(bat_pair_corr_map,callID) % only use if this callID is a key in the bat_pair_corr_map
                    % initialize values as empty
                    [l2lCorr(callID),c2lDist(callID),l2lDist(callID),c2lCorr(callID)] = deal([]);
                    
                    % get the bat IDs of the bats that produced calls in
                    % this call bout and deal with the case of overlapping
                    % calls made by multiple bats
                    calling_bat_nums = cData('callID',included_call_map(callID)').batNum;
                    if any(cellfun(@iscell,calling_bat_nums))
                        calling_bat_nums = [calling_bat_nums{:}];
                    end
                    calling_bat_nums = reshape(unique(calling_bat_nums),1,[]);
                    
                    % get this call's distance map and correlation map
                    current_call_dist_map = call_dist_by_calls(callID);
                    current_corr_map = bat_pair_corr_map(callID);
                    
                    % only use bat pairs shared between this call's
                    % distance and correlation map
                    current_bat_pairs = intersect(current_call_dist_map.keys,current_corr_map.keys);
                    % select which of these bat pairs to look at for this
                    % call
                    switch call_pair_type
                        case 'l2l' % only use bat pairs composed of bats that ~didn't~ produce this call
                            bat_pair_idx = cellfun(@(bPair) ~any(contains(bPair,calling_bat_nums)),current_bat_pairs);
                        case 'c2l' % only use bat pairs composed of bat that ~did~ produce this call
                            bat_pair_idx = cellfun(@(bPair) any(contains(bPair,calling_bat_nums)),current_bat_pairs);
                    end
                    current_bat_pairs = current_bat_pairs(bat_pair_idx);
                    
                    % iterate over the bat pairs we're using for this call
                    for batPair = current_bat_pairs
                        % get the correlation values in time around this
                        % call for this bat pair
                        batCorr = current_corr_map(batPair{1});
                        switch call_pair_type
                            case 'l2l' % we need 4 values here
                                l2lCorr(callID) = [l2lCorr(callID) nanmean(batCorr(t_idx))]; % get the average correlation between listeners during timeLims
                                l2lDist(callID) = [l2lDist(callID) current_call_dist_map(batPair{1})]; % get the average distance between listeners during timeLims
                                c2lDist = get_c2l_for_l2l(batPair,calling_bat_nums,current_call_dist_map,c2lDist); % get the average distance between the two pairs of listener-caller
                            case 'c2l' % here we only need 2
                                if any(ismember(calling_bat_nums,cData.batNums))
                                    c2lDist(callID) = [c2lDist(callID) current_call_dist_map(batPair{1})]; % get the average correlation between caller and listener during timeLims
                                    c2lCorr(callID) = [c2lCorr(callID) nanmean(batCorr(t_idx))]; % get the average distance between caller and listener during timeLims
                                end
                        end
                    end
                    k = k + 1;
                end
            end
            
            % organize into a table for linear model fitting
            switch call_pair_type
                
                case 'l2l'
                    l2lCorr_values = l2lCorr.values;l2lCorr_values = [l2lCorr_values{:}];
                    c2lDist_values = c2lDist.values;c2lDist_values = [c2lDist_values{:}];
                    l2lDist_values = l2lDist.values;l2lDist_values = [l2lDist_values{:}];
                    predTable = table(l2lDist_values',nanmean(c2lDist_values)',l2lCorr_values','VariableNames',{'l2l','c2l','corr'});
                case 'c2l'
                    c2lCorr_values = l2lCorr.values;c2lCorr_values = [c2lCorr_values{:}];
                    c2lDist_values = c2lDist.values;c2lDist_values = [c2lDist_values{:}];
                    predTable = table(c2lDist_values',c2lCorr_values','VariableNames',{'c2l','corr'});
                    
            end
            % organize into a structure for output
            call_dist_and_corr = struct('l2lCorr',l2lCorr,'c2lDist',c2lDist,'l2lDist',l2lDist,'c2lCorr',c2lCorr);
        end
        function predTable = get_session_call_dist_and_corr(pd,cData,bat_pair_corr_info,varargin)
            %% gets pairwise distance averaged across each session and average interbrain correlation across calls on the corresponding day
            % NOTE: This function is meant to analyze the relationship
            % between the average pairwise distance between bats in one
            % session and the average interbrain correlation around calls
            % in the same or another session.
            % Inputs:
            % cData: callData object for either the vocal or social session
            % bat_pair_corr_info: output of 'calculate_all_cross_brain_lfp_corr'
            % which contains the instantaneous inter-brain correlation
            % between bats around calls either during the vocal or social
            % session
            % sessionSelection: Which part of the social session to look
            % at: "all", "food", or "social"
            % excl_bat_nums: which bats to exclude
            % exclDates: which expDates to exclude
            % included_call_type: for each bat pair, which types of calls
            % to average over. Can be "all" "calling" or "listening"
            % tLim: time limits relative to call onset over which to
            % average inter-brain correlation
            % minCalls: minimum number of calls necessary to average over,
            % otherwise value is set to NaN (default 10).
            % calls_in_cluster: Flag indicating whether to restrict calls
            % (in the social session) to those that occurred at most
            % "clusterThresh" distance away
            % clusterThresh: threshold for calls occuring "in cluster"
            % Ouputs:
            % predTable: a table that can input into 'fitlm' to test the
            % relationship between distance and interbrain correlation.
            
            pnames = {'sessionSelection','excl_bat_nums','exclDates','included_call_type','tLim','minCalls','calls_in_cluster','clusterThresh'};
            dflts  = {'social','11682',datetime(2020,8,3),'all',[-0.3 0.3],10,false,100};
            [sessionSelection,excl_bat_nums,exclDates,included_call_type,tLim,minCalls,calls_in_cluster,clusterThresh] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            [bat_pair_corr_map, included_call_map] = bat_pair_corr_info_2_map(bat_pair_corr_info,pd);
            
            % if we want to restrict getting interbrain correlation values
            % to calls where bats were within cluster distance of each, get
            % the map of pairwise distance indexed by calls (only to be
            % used when analyzing calls made in the social session)
            if calls_in_cluster
                callDist = pd.get_call_dist(cData,'inter_call_int',0);
                call_dist_by_calls = pd.collect_by_calls(callDist);
            end
            
            % get list of callIDs to iterate over
            all_call_IDs = bat_pair_corr_map.keys;
            all_call_IDs = [all_call_IDs{:}];
            
            % get the list of bat pairs, we're assuming here that the list
            % doesn't change (if that's not right, then the variable
            % batCats won't be right either)
            bat_pair_keys = keys(bat_pair_corr_map(all_call_IDs(1)));
            bat_pair_keys = bat_pair_keys(~contains(bat_pair_keys,excl_bat_nums));
            
            % get the list of expDates to iterate over
            used_exp_date_strs = pd.batPos.keys;
            used_exp_dates = pd.expstr2datetime(used_exp_date_strs);
            used_exp_dates = setdiff(used_exp_dates,exclDates);
            % get the time index around calls to average over
            [~,t_idx] = inRange(bat_pair_corr_info.time,tLim);
            
            expDist = containers.Map('KeyType','char','ValueType','any');
            bat_call_corr = containers.Map('KeyType','char','ValueType','any');
            
            for expDate = used_exp_dates
                % get the expDate string used for keys
                current_date_key = pd.datetime2expstr(expDate);
                % get the pairwise distance between all bats for this day
                current_exp_dist = pd.get_pairwise_dist(expDate,'sessionSelection',sessionSelection);
                % get a (unordered) list of all callIDs made on this day
                exp_call_IDs = cData('expDay',expDate).callID';
                % only use call for which interbrain corr also exists
                exp_call_IDs = intersect(exp_call_IDs,all_call_IDs);
                % get the (ordered) list of which bats produced these calls
                exp_calling_bat_nums = cellfun(@(callID) cData('callID',callID).batNum,num2cell(exp_call_IDs));
                nCalls = length(exp_call_IDs);
                % only use the bat pairs we're interested in and that were
                % present on this day
                current_bat_pairs = intersect(bat_pair_keys,current_exp_dist.keys);
                nPair = length(current_bat_pairs);
                % get the average pairwise distance on this day for each
                % bat pair
                expDist(current_date_key) = cellfun(@(key) nanmean(current_exp_dist(key)),current_bat_pairs);
                
                current_bat_call_corr = nan(1,nPair);
                bat_pair_k = 1;
                for bPair = current_bat_pairs
                    bPair_split = strsplit(bPair{1},'-'); % split up this string to check for included_call_type
                    used_call_IDs = []; % keep track of which call bouts have already been looked at
                    call_k = 1;
                    
                    switch included_call_type
                        case 'all' % use all call bouts on this day
                            callIdx = true(1,nCalls);
                        case 'calling' % use only calls which one member of this bat pair produced
                            callIdx = cellfun(@(bNum) any(contains(bPair_split,bNum)),exp_calling_bat_nums);
                        case 'listening' % use only calls which neither member of this bat pair produced
                            callIdx = cellfun(@(bNum) ~any(contains(bPair_split,bNum)),exp_calling_bat_nums);
                    end
                    
                    % initialize list of correlation values during calls
                    % that we will average over for this pair and day
                    current_call_corr = nan(1,sum(callIdx));
                    
                    for callID = exp_call_IDs(callIdx) % iterate over the calls made on this day which are "included"
                        if calls_in_cluster % if true, restrict calls to those made when this bat pair were at most clusterThresh apart from each other
                            if isKey(call_dist_by_calls,callID)
                                current_call_dist = call_dist_by_calls(callID);
                                useCall = current_call_dist(bPair{1}) < clusterThresh;
                            else % if we don't have a distance for this call, don't use
                                useCall = false;
                            end
                        else
                            useCall = true;
                        end
                        
                        if ~ismember(callID,used_call_IDs) && useCall % use this call if we haven't looked at other calls in this bout and if it meets the cluster restrictions
                            bpCorr = bat_pair_corr_map(callID); % get this calls corr map indexed by bat pair
                            bpCorr = bpCorr(bPair{1}); % get the time varying value of the interbrain correlation for this call for this pair
                            current_call_corr(call_k) = nanmean(bpCorr(t_idx)); % average over the allotted time window
                            used_call_IDs = [used_call_IDs included_call_map(callID)]; %#ok<AGROW> % add all calls in this bout to the "used" list so we don't double count bouts
                            call_k = call_k + 1;
                        end
                    end
                    % only average over these calls if there were at least
                    % minCalls made during this day
                    if call_k - 2 > minCalls % subtract two because we 1-indexed and because of the trailing +1
                        current_bat_call_corr(bat_pair_k) = nanmean(current_call_corr);
                    end
                    bat_pair_k = bat_pair_k + 1;
                end
                % stores the list of average correlation values in the same
                % order as the average distance values to allow for
                % correlation
                bat_call_corr(current_date_key) = current_bat_call_corr;
            end
            % collect distance and correlation values into a table that can
            % be input into fitlm
            distValues = expDist.values;
            corrValues = bat_call_corr.values;
            % produce a column for the table that indicates which bat pair
            % that value came from in order to account for varying levels
            % of correlation across bat pairs in the regression (these
            % values will not be right if different groups of bats were
            % present across days)
            batCats = cellfun(@(corr) categorical(1:length(corr)),corrValues,'un',0);
            % also keep track of expDates
            all_exp_dates = cellfun(@(expDate,corr) repmat(expDate,1,length(corr)),num2cell(used_exp_dates),corrValues,'un',0);
            
            distValues = [distValues{:}];
            corrValues = [corrValues{:}];
            batCats = [batCats{:}];
            all_exp_dates = [all_exp_dates{:}];
            predTable = table(distValues',corrValues',batCats',all_exp_dates','VariableNames',{'pos','corr','bat','date'});
            
        end
        function [bat_pair_corr_map, included_call_map] = bat_pair_corr_info_2_map(bpci,pd)
            
            f_k = 3;
            nCalls = length(bpci.all_included_call_nums);
            bat_pair_keys = pd.get_pair_keys(str2double(bpci.all_bat_pairs));
            bat_pair_corr_map = containers.Map('KeyType','double','ValueType','any');
            included_call_map = containers.Map('KeyType','double','ValueType','any');
            
            for call_k = 1:nCalls
                batCorr = squeeze(bpci.bat_pair_corr(call_k,:,f_k,:));
                current_call_IDs = bpci.all_included_call_nums{call_k};
                current_call_map = containers.Map('KeyType','char','ValueType','any');
                
                bat_pair_k = 1;
                for bPair = bat_pair_keys'
                    current_call_map(bPair{1}) = batCorr(bat_pair_k,:);
                    bat_pair_k = bat_pair_k + 1;
                end
                
                for callID = current_call_IDs
                    bat_pair_corr_map(callID) = current_call_map;
                    included_call_map(callID) = current_call_IDs;
                end
            end
            
        end
        
        % functions to analyze full session lfp power and position/distance
        function [lfpResample, posResample, video_t_rs, lfp_fs] = get_aligned_lfp_pos(pd,expDate,varargin)
            pnames = {'sessionSelection'};
            dflts  = {'social'};
            [sessionSelection] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            lfpData = load_all_session_lfp_power(pd,expDate);
            [lfpResample, posResample, video_t_rs, lfp_fs] = align_lfp_pos(pd,lfpData,expDate,'sessionSelection',sessionSelection);
        end
        function [pairDist, pairCorr, video_t_rs, lfp_fs] = get_aligned_corr_dist(pd,expDate,varargin)
            
            pnames = {'lfp_fill_win','corr_smooth_scale','sessionSelection'};
            dflts  = {10,100,'social'};
            [lfp_fill_win, corr_smooth_scale, sessionSelection] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            [lfpResample, posResample, video_t_rs, lfp_fs] = get_aligned_lfp_pos(pd,expDate,'sessionSelection',sessionSelection);
            
            used_bat_nums = lfpResample.keys;
            used_bat_nums = [used_bat_nums{:}];
            % first enumerate all possible bat pairs
            batPairs = get_bat_pairs(pd,'expDate',expDate,'used_bat_nums',used_bat_nums);
            
            % get the distance between each of those pairs
            nPairs = size(batPairs,1);
            [pairDist, pairCorr] = deal(cell(1,nPairs));
            
            for bat_pair_k = 1:nPairs
                pairDist{bat_pair_k} = vecnorm(posResample(batPairs(bat_pair_k,1)) - posResample(batPairs(bat_pair_k,2)),2,2);
                
                current_lfp_data = cellfun(@(bNum) nanmedian(lfpResample(bNum),1)',num2cell(batPairs(bat_pair_k,:)),'un',0);
                current_lfp_data = fillmissing([current_lfp_data{:}],'movmean',lfp_fill_win);
                pairCorr{bat_pair_k} = movCorr(current_lfp_data(:,1),current_lfp_data(:,2),corr_smooth_scale,0);
            end
            
            % convert the array of bat pairs into a list of strings of the
            % form '[batNum1]-[batNum2]' to use as keys into the Map of
            % pairwise distances
            bat_pair_keys = pd.get_pair_keys(batPairs);
            pairDist = containers.Map(bat_pair_keys,pairDist);
            pairCorr = containers.Map(bat_pair_keys,pairCorr);
            
        end
        function [all_pair_dist, all_pair_corr] = get_all_session_corr_dist(pd,varargin)
            
            pnames = {'excl_bat_nums','corr_smooth_scale','exclDates'};
            dflts  = {'11682',100,{'08052020','08032020'}};
            [excl_bat_nums, corr_smooth_scale, exclDates] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            batNums = str2double(setdiff(pd.batNums,excl_bat_nums));
            batPairs = pd.get_bat_pairs('used_bat_nums',batNums);
            bat_pair_keys = pd.get_pair_keys(batPairs);
            
            all_pair_corr = containers.Map('KeyType','char','ValueType','any');
            all_pair_dist = containers.Map('KeyType','char','ValueType','any');
            
            expDates = pd.batPos.keys;
            expDates = setdiff(expDates,exclDates);
            expDates = pd.expstr2datetime(expDates);
            
            for expDate = expDates
                [pairDist, pairCorr] = get_aligned_corr_dist(pd,expDate,'corr_smooth_scale',corr_smooth_scale);
                for batPair = bat_pair_keys'
                    if isKey(pairCorr,batPair{1})
                        if ~isKey(all_pair_corr,batPair{1})
                            all_pair_corr(batPair{1}) = pairCorr(batPair{1});
                            all_pair_dist(batPair{1}) = pairDist(batPair{1});
                        else
                            all_pair_corr(batPair{1}) = [all_pair_corr(batPair{1}); pairCorr(batPair{1})];
                            all_pair_dist(batPair{1}) = [all_pair_dist(batPair{1}); pairDist(batPair{1})];
                        end
                    end
                end
            end
        end
    end
end

function c2lDist = get_c2l_for_l2l(batPair,calling_bat_nums,current_call_dist_map,c2lDist)

listening_bat_nums = strsplit(batPair{1},'-');
listening_bat_dist = nan(2,1);
listen_bat_k = 1;
for listenBat = listening_bat_nums
    current_listen_bat_dist = nan(1,length(calling_bat_nums));
    call_bat_k = 1;
    for callBat = calling_bat_nums
        current_bat_pair_key = strjoin(sort([callBat listenBat]),'-');
        current_listen_bat_dist(call_bat_k) = current_call_dist_map(current_bat_pair_key);
        call_bat_k = call_bat_k + 1;
    end
    listening_bat_dist(listen_bat_k) = nanmean(current_listen_bat_dist);
    listen_bat_k = listen_bat_k + 1;
end
c2lDist(callID) = [c2lDist(callID) listening_bat_dist];
end

function pred_centroids = get_pred_centroids(LEDTracks)

%% Inputs:
% LEDTracks: results of LED tracking on an entire session (includes
% centroidLocs, predColors, and predPosterior).
% color_pred_model: classification model used to distinguish colors, used
% here as an index into the prediction posterior matrix that the model
% produces.
%% Outputs:
% pred_centroids: a [num. of frames X 2 X num. of colors] matrix of
% prediction centroid locations in X and Y for each frame and each color.

color_names = LEDTracks.color_pred_model.ClassificationSVM.ClassNames;
nColor = length(color_names);

nFrames = length(LEDTracks.centroidLocs);
pred_centroids = nan(nFrames,2,nColor);
for frame_k = 1:nFrames
    for color_k = 1:nColor
        current_color_idx = strcmp(LEDTracks.predColors{frame_k},color_names{color_k}); % which color relative to the list used in the prediction model are we looking at?
        if sum(current_color_idx) > 1 % if the model predicted more than 1 of the same color, decide which to use
            current_pred_posteriors = LEDTracks.predPosterior{frame_k}(current_color_idx,color_k);
            if length(unique(current_pred_posteriors)) == 1 % if multiple locations have the same posterior (e.g. both at maximum posterior), we can't decide between them, discard
                pred_centroids(frame_k,:,color_k) = NaN;
                continue
            else
                [~,pred_posterior_idx] = max(current_pred_posteriors); % take the prediction with the higher posterior as the correct prediction
                current_color_idx = find(current_color_idx);
                current_color_idx = current_color_idx(pred_posterior_idx);
            end
        elseif sum(current_color_idx) == 0 % if this color isn't present in the predictin mark as NaN
            pred_centroids(frame_k,:,color_k) = NaN;
            continue
        end
        pred_centroids(frame_k,:,color_k) = LEDTracks.centroidLocs{frame_k}(current_color_idx,:);
    end
end

end

function pred_centroids_clean = clean_pred_centroids(pred_centroids, varargin)
% First-pass analysis of location data taken from the raw results of the LED_tracking model for a given day.

% first, we remove douplicaitons of color per frame (if you already have the variable locs saved you can load it and it will skip this first part*)
% We then transpose the xy values to be able to empose the XY limts of cage top
% next, we filter and fill in gaps in data.
% INPUTS:
% 1)file: the dir location and and name of prediction results to load (could be locs alone)
% optional:
% 2)filtRank: the rank of the median filter used to smooth data (scalar, #frames)
% 3)movWindow: the size of the moving window for filling the gaps if useing defult (movMedian), (scalar, #frames)
% 4)GapTh: if we want to limit to only gaps of a certain size or under we specify this (scalar, #frames)
% 5)GapMethod: if we gave Gapth, we can now choose the method, e.g change to 'spline', (string)

% * the variable should be named locs in the workspace, it is a variable with 3D. the 1D(rows) is number of frames for the whole day
% the 2D(colomns) is number of colors(e.g 8), 3D is x,y values for each frame. e.g for dimentions: 144314x8x2.

pnames = {'filtRank','movWindow','fillGaps','GapMethod','remove_out_of_bounds'};
dflts  = {5,5,true,'linear_with_th',true};
[filtRank, movWindow, fillGaps, GapMethod, remove_out_of_bounds] = internal.stats.parseArgs(pnames,dflts,varargin{:});

nColor = size(pred_centroids,3);
pred_centroids_clean = nan(size(pred_centroids));

%% (2) here we rotate so we can limit x,y to max limits, med_filter the data, fill gaps.
ROI_rot = load('ROI_rot.mat'); % this struc contains the rotation matrix, limts of x and y and center of cage values
for color_k = 1:nColor % we run this analysis bat by bat (e.i color by color)
    
    % "rotate" pixels about a center point
    XYrot = ((pred_centroids(:,:,color_k) - ROI_rot.c)*ROI_rot.R') + ROI_rot.c;
    
    % enforce x,y min/max limits for stepping out of cage-top bounds
    
    
    if ~remove_out_of_bounds
        x = XYrot(:,1);
        y = XYrot(:,2);
        
        x(x<ROI_rot.xlims(1)) = ROI_rot.xlims(1);
        x(x>ROI_rot.xlims(2)) = ROI_rot.xlims(2);
        
        y(y<ROI_rot.ylims(1)) = ROI_rot.ylims(1);
        y(y>ROI_rot.ylims(2)) = ROI_rot.ylims(2);
        
        XYrot = [x y];
    else
        OOBIdx = XYrot(:,1) < ROI_rot.xlims(1) | XYrot(:,1) > ROI_rot.xlims(2) |...
            XYrot(:,2) < ROI_rot.ylims(1) | XYrot(:,2) > ROI_rot.ylims(2);
        XYrot(OOBIdx,:) = NaN;
    end
    
    XYrot = XYrot - [ROI_rot.xlims(1) ROI_rot.ylims(1)];
    
    % median filter so that the median of any segment containing NaNs is
    % the median of the non-NaN values and so that segments around the
    % edges are truncated (instead of assuming zero padding)
    XYfilt = medfilt1(XYrot,filtRank,'omitnan','truncate');
    % gap filling:
    if fillGaps
        switch GapMethod
            case 'linear_with_th' % linearly interpolate gaps shorter than GapTh samples
                gap_bw_labels = bwlabel(isnan(XYfilt(:,1))); % fast way to tag gaps by using bwlable to find "connceted" regions of logical mtraix in 1D, and lable them
                gapLabels = setdiff(unique(gap_bw_labels),0); % get a list of the uniquely labeled gaps
                gapSizes = histcounts(gap_bw_labels,[gapLabels; gapLabels(end)+1]); % count the number of samples in each gap
                
                used_gap_labels = gapLabels(gapSizes <= movWindow); % find gap segments shorter than GapTh
                missingIdx = ismember(gap_bw_labels,used_gap_labels); % get an index of samples to be filled
                XYfilt_filled = fillmissing(XYfilt,'linear','MissingLocations',missingIdx);
            case 'movmedian' % fill gaps across the whole trace using the median moving window
                XYfilt_filled = fillmissing(XYfilt,'movmedian',movWindow);
        end
    else
        XYfilt_filled = XYfilt;
    end
    
    % collecting the results:
    pred_centroids_clean(:,:,color_k) = XYfilt_filled;
end


end

function [idx_tracks,idx_ts] = align_tracks_with_frame_ts(LEDTracks,frame_ts_info)
% gets indices of timestamps and LED tracks that share the same video file
% and frame number
frame_and_file_table_data = array2table([LEDTracks.file_frame_number;LEDTracks.fileIdx]');
frame_and_file_table_ts = array2table([frame_ts_info.file_frame_number;frame_ts_info.fileIdx]');
[~,idx_tracks,idx_ts] = intersect(frame_and_file_table_data,frame_and_file_table_ts,'stable');
end

function reordered_bat_pos = reorder_bat_pos(pd,batPos,colorStrs,expDate)
% takes a nSample X 2 X nBat array of XY locations orderd according to the
% list of colors in the color prediction model, gets the batNums
% correspoding to those colors for this expDate and reorders it according
% to the order of bats present in pd.all_bat_nums

bat_color_table = get_bat_color_table(pd,expDate);

nBats = length(pd.all_bat_nums);
nSample = size(batPos,1);
reordered_bat_pos = nan(nSample,2,nBats);

for color_k = 1:length(colorStrs)
    color_bat_num = bat_color_table.batNum(strcmp(colorStrs{color_k},bat_color_table.color));
    current_bat_idx = pd.all_bat_nums == color_bat_num;
    reordered_bat_pos(:,:,current_bat_idx) = batPos(:,:,color_k);
end

end

function bat_color_table = get_bat_color_table(pd,expDate)
% returns a table of batNums and their corresponding audioLogger color
% strings used this expDate
if ~isdatetime(expDate)
    expDate = pd.expstr2datetime(expDate);
end
bat_idx = contains(pd.recLogs.Properties.VariableNames,'Bat_');
color_idx = contains(pd.recLogs.Properties.VariableNames,'Color_');
dateIdx = pd.recLogs.Date == expDate & strcmp(pd.recLogs.Session,pd.sessionType);
assert(sum(dateIdx) == 1)
rec_logs_exp = pd.recLogs(dateIdx,:);
bat_color_table = table(rec_logs_exp{1,bat_idx}',rec_logs_exp{1,color_idx}','VariableNames',{'batNum','color'});
end

function sync_bat_num = get_sync_bat_num(pd,expDate)
%% gets the bat that was chose as the 'sync' bat, which has the relevant session strings logged in its eventlog file
dateIdx = pd.recLogs.Date == expDate & strcmp(pd.recLogs.Session,pd.sessionType);

sync_logger_num = num2str(pd.recLogs.Sync_logger_num(dateIdx));

NL_idx = find(contains(pd.recLogs.Properties.VariableNames,'NL_'));
NLStrs = strsplit(pd.recLogs.Properties.VariableNames{NL_idx(pd.recLogs{dateIdx,NL_idx} == str2double(sync_logger_num))},'_');
batStr = strjoin({'Bat',NLStrs{2}},'_');
sync_bat_num = num2str(pd.recLogs.(batStr)(dateIdx));

end

function foodTime = get_food_time(pd,expDate,varargin)
%% gets the time of food delivery by looking in the eventlog file for a food related string 

pnames = {'foodStr'};
dflts  = {'banana'};
[foodStr] = internal.stats.parseArgs(pnames,dflts,varargin{:});

call_data_dir = fullfile(pd.baseDirs{1},'call_data');
event_file_dir = fullfile(pd.baseDirs{1},'event_file_data');
exp_date_str = datestr(expDate,'yyyymmdd');
% get this day's recording logs and check which session came first, social
% or vocal 
exp_rec_logs = pd.recLogs(pd.recLogs.Date == expDate & ~strcmp(pd.recLogs.Session,'playback'),:);
first_session_type = exp_rec_logs.Session{1};

% get the audio2nlg file for whichever session came first, as the first TTL
% pulse from that session is defined as time zero
switch first_session_type
    case 'vocal'
        audio2nlg_fname = fullfile(call_data_dir,[exp_date_str '_audio2nlg_fit.mat']);
    case 'social'
        audio2nlg_fname = fullfile(call_data_dir,[exp_date_str '_audio2nlg_fit_social.mat']);
end

audio2nlg = load(audio2nlg_fname);

% sync bat is he one that will have the foodstr in its log
sync_bat_num = get_sync_bat_num(pd,expDate);
event_file_fname = fullfile(event_file_dir,[sync_bat_num '_' exp_date_str '_EVENTS.mat']);
eventData = load(event_file_fname);

% get the food time, and if the foodStr couldn't be found, return NaN
foodIdx = contains(eventData.event_types_and_details,foodStr);
if any(foodIdx)
    foodTime = 1e-3*eventData.event_timestamps_usec(foodIdx) - audio2nlg.first_nlg_pulse_time;
else
    foodTime = NaN;
end

end

function [xSub,idx] = inRange(x,bounds)
% function to get values and index of  vector within a given bounds
bounds = sort(bounds);
idx = x>bounds(1) & x<= bounds(2);
xSub = x(idx);
end

function all_session_lfp_power = load_all_session_lfp_power(pd,expDate,varargin)
%% loads the pre-calculated lfp power for all bats present on a given date
% Inputs: 
% expDate: which expDate to get lfp data for
% max_artifacts_fract: maximum fraction of samples in a window in which lfp
% power was calculated to use, any window over that will be set to NaN
% freq_k: which frequency band to use

pnames = {'max_artifact_fract','freq_k'};
dflts  = {0.01,2};
[max_artifact_fract, freq_k] = internal.stats.parseArgs(pnames,dflts,varargin{:});

lfp_data_dir = fullfile(pd.baseDirs{1},'lfp_data');
if isdatetime(expDate)
    expDate = datestr(expDate,'yyyymmdd');
end

lfpFnames = dir(fullfile(lfp_data_dir,['*' expDate '_all_session_lfp_results.mat']));
for bat_k = 1:length(lfpFnames)
    % load the relevant data
    current_lfp_power = load(fullfile(lfpFnames(bat_k).folder,lfpFnames(bat_k).name),'lfpPower','winSize','lfp_power_timestamps','batNum','n_artifact_times'); 
    % remove windows that had too many artifacts detected
    current_lfp_power.lfpData = get_artifact_removed_full_session_LFP(current_lfp_power,max_artifact_fract,freq_k); 
    current_lfp_power.batNum = str2double(current_lfp_power.batNum); 
    all_session_lfp_power(bat_k) = current_lfp_power; %#ok<AGROW>
end

end

function [lfpResample, posResample, video_t_rs, lfp_fs] = align_lfp_pos(pd,lfpData,expDate,varargin)

pnames = {'sessionSelection'};
dflts  = {'social'};
[sessionSelection] = internal.stats.parseArgs(pnames,dflts,varargin{:});

exp_date_str = pd.datetime2expstr(expDate);
fTime = 1e-3*pd.foodTime(exp_date_str);
video_t = 1e-3*pd.posTS(exp_date_str);

if ~isnan(fTime)
    switch sessionSelection
        case 'social'
            session_t_idx = video_t < fTime;
        case 'food'
            session_t_idx = video_t > fTime;
        case 'all'
            session_t_idx = true(size(video_t));
    end
else
    session_t_idx = true(size(video_t));
end

roundingFactor = 100;

lfp_fs = 1/median(diff(lfpData(1).lfp_power_timestamps));
lfp_fs = round(roundingFactor*lfp_fs)/roundingFactor;

[N,D] = rat(lfp_fs/pd.video_fs);

tRange = [max(cellfun(@min,{lfpData.lfp_power_timestamps})) min(cellfun(@max,{lfpData.lfp_power_timestamps}))];

[~,lfp_t_idx] = inRange(video_t,tRange);
t_range_idx = lfp_t_idx & session_t_idx;
video_t = video_t(t_range_idx);
video_t_rs = resample(video_t,N,D,0);

pos = pd.get_pos(exp_date_str);
batNums = pos.keys;
posResample = containers.Map('KeyType','double','ValueType','any');

for bNum = batNums
    current_bat_pos = pos(bNum{1});
    posResample(bNum{1}) = resample(current_bat_pos(t_range_idx,:),N,D);
end

nBat = length(lfpData);
lfpResample = containers.Map('KeyType','double','ValueType','any');
for bat_k = 1:nBat
    current_lfp_power = lfpData(bat_k).lfpData;
    nChannel = size(current_lfp_power,1);
    current_lfp_resample = nan(nChannel,length(video_t_rs));
    for ch_k = 1:nChannel
        current_lfp_resample(ch_k,:) = interp1(lfpData(bat_k).lfp_power_timestamps,current_lfp_power(ch_k,:),video_t_rs,'linear');
    end
    lfpResample(lfpData(bat_k).batNum) = current_lfp_resample;
end

video_t_rs = 1e3*video_t_rs;

end

function lfpPower_artifact_removed = get_artifact_removed_full_session_LFP(lfpData,max_artifact_frac,freq_k)

artifact_chunks = lfpData.n_artifact_times/lfpData.winSize > max_artifact_frac;
lfpPower_artifact_removed = lfpData.lfpPower(:,:,freq_k);
lfpPower_artifact_removed(artifact_chunks) = NaN;

end