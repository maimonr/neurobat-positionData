classdef positionData < ephysData
    properties
        gap_fill_window_s = 2
        video_fs = 20
        all_bat_nums = [11636,11682,13688,14612,14654,14798,60064,71216,71335]
        sessionType = 'social'
        callOffset = 2
        
        batPos
        posTS
        tracking_data_dir
        video_data_dir
        call_data_dir
        
    end
    
    methods
        
        function pd = positionData(varargin)
            
            pnames = {'used_exp_dates'};
            dflts  = {[]};
            [used_exp_dates] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            pd = pd@ephysData('adult_social');
            pd.tracking_data_dir = fullfile(pd.baseDirs{1},'tracking_data');
            pd.video_data_dir = fullfile(pd.baseDirs{1},'video_data');
            pd.call_data_dir = fullfile(pd.baseDirs{1},'call_data');
            
            mov_window_samples = pd.gap_fill_window_s*pd.video_fs;
            
            led_track_fnames = dir(fullfile(pd.tracking_data_dir,['LEDtracking_pred_' pd.sessionType '*.mat']));
            fnameSplit = arrayfun(@(x) strsplit(x.name,'_'),led_track_fnames,'un',0);
            exp_date_strs = cellfun(@(x) x{end}(1:end-4),fnameSplit,'un',0);
            expDates = pd.expstr2datetime(exp_date_strs);
            
            if ~isempty(used_exp_dates)
                [~,used_date_idx] = ismember(expDates,used_exp_dates);
                expDates = expDates(used_date_idx);
                exp_date_strs = exp_date_strs(used_date_idx);
                led_track_fnames = led_track_fnames(used_date_idx);
            end
            
            nExp = length(exp_date_strs);
            used_exp_idx = true(1,nExp);
            [bat_pos_cell, pos_ts_cell] = deal(cell(1,nExp));
            
            for exp_k = 1:nExp
                led_track_fname = fullfile(led_track_fnames(exp_k).folder,led_track_fnames(exp_k).name);
                exp_date_str = exp_date_strs{exp_k};
                frame_ts_info_fname = fullfile(pd.video_data_dir,[exp_date_str '_color_frame_timestamps_info_social.mat']);
                
                if ~isfile(frame_ts_info_fname)
                    used_exp_idx(exp_k) = false;
                    continue
                end
                s = load(frame_ts_info_fname);
                frame_ts_info = s.frame_ts_info;
                LEDTracks = load(led_track_fname);
                
                [idx_tracks,idx_ts] = align_tracks_with_frame_ts(LEDTracks,frame_ts_info);
                timestampsNlg = frame_ts_info.timestamps_nlg(idx_ts);
                
                assert(round(1e3/median(diff(timestampsNlg))) == pd.video_fs)
                
                predCentroids = get_pred_centroids(LEDTracks);
                predCentroids = predCentroids(idx_tracks,:,:);
                predCentroids = clean_pred_centroids(predCentroids,'GapMethod','movmedian','movWindow',mov_window_samples);
                
                colorStrs = LEDTracks.color_pred_model.ClassificationSVM.ClassNames;
                predCentroids = reorder_bat_pos(pd,predCentroids,colorStrs,expDates(exp_k));
                predCentroids = squeeze(num2cell(predCentroids,[1 2]));
                
                bat_pos_cell{exp_k} = containers.Map(pd.all_bat_nums,predCentroids);
                pos_ts_cell{exp_k} = timestampsNlg';
                
            end
            
            pd.batPos = containers.Map(exp_date_strs(used_exp_idx),bat_pos_cell(used_exp_idx));
            pd.posTS = containers.Map(exp_date_strs(used_exp_idx),pos_ts_cell(used_exp_idx));
            
        end
        function pos = get_pos(pd,expDate,batNum)
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
            if isdatetime(expDate)
                expDate = datestr(expDate,'mmddyyyy');
            end
            t = pd.posTS(expDate);
        end
        function dist = get_dist(pd,batNums,expDate)
            dist = vecnorm((get_pos(pd,expDate,batNums(1)) - get_pos(pd,expDate,batNums(2))),2,2);
        end
        function pairDist = get_pairwise_dist(pd,expDate,varargin)
            pnames = {'used_bat_nums'};
            dflts  = {[]};
            [used_bat_nums] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            batPairs = get_bat_pairs(pd,'expDate',expDate,'used_bat_nums',used_bat_nums);
            nPairs = size(batPairs,1);
            pairDist = cell(1,nPairs);
            for bat_pair_k = 1:nPairs
                pairDist{bat_pair_k} = get_dist(pd,batPairs(bat_pair_k,:),expDate);
            end
            bat_pair_keys = pd.get_pair_keys(batPairs);
            pairDist = containers.Map(bat_pair_keys,pairDist);
            
        end
        function batPairs = get_bat_pairs(pd,varargin)
            
            pnames = {'expDate','used_bat_nums'};
            dflts  = {[],[]};
            [expDate,bat_num_list] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            used_bat_nums = pd.all_bat_nums;
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
            bat_pair_keys = cellfun(@(batPair) num2str(batPair,'%d-%d'),num2cell(batPairs,2),'UniformOutput',false);
        end
        function batPairs = get_key_pairs(~,bat_pair_keys)
            batPairs = cellfun(@(pairKey) str2double(strsplit(pairKey,'-')),bat_pair_keys,'un',0);
            batPairs = vertcat(batPairs{:});
        end
        function callPos = get_call_pos(pd,cData,varargin)
            
            pnames = {'inter_call_int','expDates'};
            dflts  = {[],[]};
            [inter_call_int,expDates] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            if isempty(expDates)
                exp_date_strs = pd.batPos.keys;
                expDates = pd.expstr2datetime(exp_date_strs);
            else
                exp_date_strs = pd.datetime2expstr(expDates);
            end
            nExp = length(expDates);
            callPos = cell(1,nExp);
            for exp_k = 1:nExp
                pos = pd.get_pos(expDates(exp_k));
                callTimes = cData('expDay',expDates(exp_k)).callPos;
                callTimes = 1e3*callTimes(:,1);
                calling_bat_nums = cData('expDay',expDates(exp_k)).batNum';
                
                if ~isempty(inter_call_int)
                    callIdx = [Inf; diff(callTimes)] > inter_call_int;
                    callTimes = callTimes(callIdx);
                    calling_bat_nums = calling_bat_nums(callIdx);
                end
                
                t = pd.get_time(expDates(exp_k));
                call_time_offset = 1e3.*pd.callOffset.*[-1 1];
                nCall = length(callTimes);
                nBat = length(pd.all_bat_nums);
                call_pos_idx = cell(1,nCall);
                current_call_pos = repmat({nan(nCall,2)},1,nBat);
                
                for call_k = 1:nCall
                    current_call_time = callTimes(call_k) + call_time_offset;
                    [~,call_pos_idx{call_k}] = inRange(t, current_call_time);
                end
                
                for bat_k = 1:nBat
                    current_bat_pos = pos(pd.all_bat_nums(bat_k));
                    current_bat_pos = cellfun(@(idx) nanmean(current_bat_pos(idx,:)),call_pos_idx,'un',0);
                    current_call_pos{bat_k} = struct('pos',current_bat_pos,'caller',calling_bat_nums);
                end
                
                callPos{exp_k} = containers.Map(pd.all_bat_nums,current_call_pos);
            end
            callPos = containers.Map(exp_date_strs,callPos);
        end
        function callDist = get_call_dist(pd,cData,varargin)
            
            pnames = {'inter_call_int','expDates'};
            dflts  = {[],[]};
            [inter_call_int,expDates] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            if isempty(expDates)
                exp_date_strs = pd.batPos.keys;
                expDates = pd.expstr2datetime(exp_date_strs);
            else
                exp_date_strs = pd.datetime2expstr(expDates);
            end
            
            call_bat_pos = get_call_pos(pd,cData,'expDate',expDates,'inter_call_int',inter_call_int);
            
            nExp = length(expDates);
            callDist = cell(1,nExp);
            for exp_k = 1:nExp
                batPairs = get_bat_pairs(pd,'expDate',expDates(exp_k));
                nPairs = size(batPairs,1);
                current_bat_pos = call_bat_pos(exp_date_strs{exp_k});
                bat_call_dist = cell(1,nPairs);
                for bat_pair_k = 1:nPairs
                    current_call_bat_pos = cellfun(@(bat) current_bat_pos(bat),num2cell(batPairs(bat_pair_k,:)),'un',0);
                    calling_bat_nums = cellfun(@(bat) {bat.caller},current_call_bat_pos,'un',0);
                    assert(all(~cellfun(@ischar,calling_bat_nums{1}) | strcmp(calling_bat_nums{:})))
                    
                    current_call_bat_pos = cellfun(@(bat) vertcat(bat.pos),current_call_bat_pos,'un',0);
                    current_call_dist = vecnorm(current_call_bat_pos{1} - current_call_bat_pos{2},2,2)';
                    bat_call_dist{bat_pair_k} = struct('pos',num2cell(current_call_dist),'caller',calling_bat_nums{1});
                end
                bat_pair_keys = pd.get_pair_keys(batPairs);
                callDist{exp_k} = containers.Map(bat_pair_keys,bat_call_dist);
            end
            callDist = containers.Map(exp_date_strs,callDist);
        end
        
        function exp_date_strs = datetime2expstr(~,expDates)
            exp_date_strs = cellfun(@(expDate) datestr(expDate,'mmddyyyy'),num2cell(expDates),'un',0);
        end
        function expDates = expstr2datetime(~,exp_date_strs)
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

frame_and_file_table_data = array2table([LEDTracks.file_frame_number;LEDTracks.fileIdx]');
frame_and_file_table_ts = array2table([frame_ts_info.file_frame_number;frame_ts_info.fileIdx]');
[~,idx_tracks,idx_ts] = intersect(frame_and_file_table_data,frame_and_file_table_ts,'stable');


end

function reordered_bat_pos = reorder_bat_pos(pd,batPos,colorStrs,expDate)

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
bat_idx = contains(pd.recLogs.Properties.VariableNames,'Bat_');
color_idx = contains(pd.recLogs.Properties.VariableNames,'Color_');
dateIdx = pd.recLogs.Date == expDate & strcmp(pd.recLogs.Session,pd.sessionType);
assert(sum(dateIdx) == 1)
rec_logs_exp = pd.recLogs(dateIdx,:);
bat_color_table = table(rec_logs_exp{1,bat_idx}',rec_logs_exp{1,color_idx}','VariableNames',{'batNum','color'});
end

function [xSub,idx] = inRange(x,bounds)
bounds = sort(bounds);
idx = x>bounds(1) & x<= bounds(2);
xSub = x(idx);
end



