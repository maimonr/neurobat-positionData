function clustIdx = fill_gaps_clusters(clustIdx,batPos)

gapWin = 200;
distThresh = 15;

missingFrames = all(isnan(clustIdx),2);

gapOn = find(diff(missingFrames) == 1);
gapOff = find(diff(missingFrames) == -1)+1;

assert(abs(length(gapOn) - length(gapOff)) <= 1);

if length(gapOn) > length(gapOff)
    gapOn = gapOn(1:end-1);
elseif length(gapOff) > length(gapOn)
    gapOff = gapOff(2:end);
end

nGap = length(gapOn);
for gap_k = 1:nGap
    
    pre_gap_idx = gapOn(gap_k)-gapWin:gapOn(gap_k);
    post_gap_idx = gapOff(gap_k):gapOff(gap_k)+gapWin;
    
    [pre_gap_clustPos,pre_gap_clustNums] = get_clust_pos(clustIdx(pre_gap_idx,:),batPos(:,:,pre_gap_idx));
    [post_gap_clustPos,post_gap_clustNums] = get_clust_pos(clustIdx(post_gap_idx,:),batPos(:,:,post_gap_idx));
    
    for clust_k = 1:length(pre_gap_clustNums)
        d = vecnorm(pre_gap_clustPos(clust_k,:) - post_gap_clustPos,2,2);
        if sum(d < distThresh) == 1
           current_cluster = pre_gap_clustNums(clust_k);
           cluster_to_replace = post_gap_clustNums(d < distThresh);
           clustIdx(clustIdx == cluster_to_replace) = current_cluster;
        end
    end
    
end

end

function [clustPos,clustNums] = get_clust_pos(clustIdx_gap,batPos_gap)

clustNums = unique(clustIdx_gap(:));
clustNums = clustNums(clustNums > 0);

clustPos = nan(length(clustNums),2);

for clust_k = 1:length(clustNums)
    clustNum = clustNums(clust_k);
    batIdx = any(clustIdx_gap == clustNum,1);
    clustPos(clust_k,:) = squeeze(nanmedian(batPos_gap(batIdx,:,:),[1 3]));
end

end
