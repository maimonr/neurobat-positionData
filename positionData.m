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
        
        batPos % XY positions of each bat in a nested map, indexed by expDate and then batNum
        posTS % timestamps of positions in a map indexed by expDate
        tracking_data_dir
        video_data_dir
        call_data_dir
        
    end
    
    methods
        
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
            pd.call_data_dir = fullfile(pd.baseDirs{1},'call_data');
            
            mov_window_samples = pd.gap_fill_window_s*pd.video_fs; % get number of frames over which to perform moving window cleaning of positions
            
            % get a list of all LEDtracking results and the corresponding
            % expDates
            led_track_fnames = dir(fullfile(pd.tracking_data_dir,['LEDtracking_pred_' pd.sessionType '*.mat'])); 
            fnameSplit = arrayfun(@(x) strsplit(x.name,'_'),led_track_fnames,'un',0);
            exp_date_strs = cellfun(@(x) x{end}(1:end-4),fnameSplit,'un',0);
            expDates = pd.expstr2datetime(exp_date_strs);
            
            if ~isempty(used_exp_dates) % if supplied, limit expDates to those supplied
                [~,used_date_idx] = ismember(expDates,used_exp_dates);
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
                
                % load frame_ts_info and LEDtracks
                s = load(frame_ts_info_fname);
                frame_ts_info = s.frame_ts_info;
                LEDTracks = load(led_track_fname);
                
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
            
        end
        function pos = get_pos(pd,expDate,batNum)
            % utility function to get bat position(s) for a given expDate
            % and, optionally, for a given batNum. expDate can be either a
            % datetime or a date string. Without a batNum, this returns
            % the map of bat positions for the expDate, with a batNum this
            % returns the array of positions for that bat on that date
            if isdatetime(expDate)
                expDate = datestr(expDate,'mmddyyyy');
            end
            pos = pd.batPos(expDate);
            if nargin > 2
                if ~isnumeric(batNum)
                    batNum = str2double(batNum);
                end
                pos = pos(batNum);
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
        function dist = get_dist(pd,batNums,expDate)
            % gets the frame by frame euclidean distance between two bats
            % on a given expDate
            % Inputs: 
            % batNums: 1 x 2 array of batNums
            % expDate: datetime or date string
            dist = vecnorm((get_pos(pd,expDate,batNums(1)) - get_pos(pd,expDate,batNums(2))),2,2);
        end
        function pairDist = get_pairwise_dist(pd,expDate,varargin)
            % gets pairwise distance between all pairs of bats. 
            % Inputs:
            % used_bat_nums: optional list to limit which bats to look at
            % Outputs:
            % pairDist: a new map indexed by strings of the form
            % '[batNum1]-[batNum2]' whose values are the frame by frame
            % distance between that pair.
            
            pnames = {'used_bat_nums'};
            dflts  = {[]};
            [used_bat_nums] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            % first enumerate all possible bat pairs 
            batPairs = get_bat_pairs(pd,'expDate',expDate,'used_bat_nums',used_bat_nums);
            
            % get the distance between each of those pairs
            nPairs = size(batPairs,1);
            pairDist = cell(1,nPairs);
            for bat_pair_k = 1:nPairs
                pairDist{bat_pair_k} = get_dist(pd,batPairs(bat_pair_k,:),expDate);
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
            bat_pair_keys = cellfun(@(batPair) num2str(batPair,'%d-%d'),num2cell(batPairs,2),'UniformOutput',false);
        end
        function batPairs = get_key_pairs(~,bat_pair_keys)
            % utility function to convert a list of bat pair strings back
            % to a nPair x 2 array of batNums
            batPairs = cellfun(@(pairKey) str2double(strsplit(pairKey,'-')),bat_pair_keys,'un',0);
            batPairs = vertcat(batPairs{:});
        end
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
            
            pnames = {'inter_call_int','expDates'};
            dflts  = {[],[]};
            [inter_call_int,expDates] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
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
                
                % if provided, limit to calls separated by minimum inter
                % call interval
                if ~isempty(inter_call_int)
                    callIdx = [Inf; diff(callTimes)] > inter_call_int;
                    callTimes = callTimes(callIdx);
                    calling_bat_nums = calling_bat_nums(callIdx);
                end
                
                % get the frame timestamps for this expDate
                t = pd.get_time(expDates(exp_k));
                % set up +/- offset in ms
                call_time_offset = 1e3.*pd.callOffset.*[-1 1];
                
                % here we'll iterate over all calls and then over bats
                nCall = length(callTimes);
                nBat = length(pd.all_bat_nums);
                call_pos_idx = cell(1,nCall);
                current_call_pos = repmat({nan(nCall,2)},1,nBat);
                
                % for each call get the corresponding index into the
                % position data over which to average
                for call_k = 1:nCall
                    current_call_time = callTimes(call_k) + call_time_offset;
                    [~,call_pos_idx{call_k}] = inRange(t, current_call_time);
                end
                
                % for each bat, get its average position for all calls and
                % store as a struct with fields 'pos' and 'caller' which
                % contain the given bat's position and which bat made the
                % call (we might want to save also the unique callID here)
                for bat_k = 1:nBat
                    current_bat_pos = pos(pd.all_bat_nums(bat_k));
                    current_bat_pos = cellfun(@(idx) nanmean(current_bat_pos(idx,:)),call_pos_idx,'un',0);
                    current_call_pos{bat_k} = struct('pos',current_bat_pos,'caller',calling_bat_nums);
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
            pnames = {'inter_call_int','expDates'};
            dflts  = {[],[]};
            [inter_call_int,expDates] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            % get expDates in datetime and corresponding date strings
            if isempty(expDates)
                exp_date_strs = pd.batPos.keys;
                expDates = pd.expstr2datetime(exp_date_strs);
            else
                exp_date_strs = pd.datetime2expstr(expDates);
            end
            
            % get bat positions around calls across expDates
            call_bat_pos = get_call_pos(pd,cData,'expDate',expDates,'inter_call_int',inter_call_int);
            
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
                    calling_bat_nums = cellfun(@(bat) {bat.caller},current_call_bat_pos,'un',0);
                    assert(all(~cellfun(@ischar,calling_bat_nums{1}) | strcmp(calling_bat_nums{:})))
                    
                    % get the array of positions for both bat in this pair
                    current_call_bat_pos = cellfun(@(bat) vertcat(bat.pos),current_call_bat_pos,'un',0);
                    % calculate the distance between this pair of bats
                    current_call_dist = vecnorm(current_call_bat_pos{1} - current_call_bat_pos{2},2,2)';
                    % save as a struct with fields 'dist' and 'caller'
                    bat_call_dist{bat_pair_k} = struct('dist',num2cell(current_call_dist),'caller',calling_bat_nums{1});
                end
                % create map of structs indexed by bat pair strings
                bat_pair_keys = pd.get_pair_keys(batPairs);
                callDist{exp_k} = containers.Map(bat_pair_keys,bat_call_dist);
            end
            % create a map of maps indexed by expDate
            callDist = containers.Map(exp_date_strs,callDist);
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
        
    end
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

function [xSub,idx] = inRange(x,bounds)
% function to get values and index of  vector within a given bounds
bounds = sort(bounds);
idx = x>bounds(1) & x<= bounds(2);
xSub = x(idx);
end



