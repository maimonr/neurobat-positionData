function clustIdx = track_clusters(clustIdx)
clustIdx_orig = clustIdx;
dead_clusters = [];
new_clust_num = max(clustIdx,[],'all') + 1;
nT = size(clustIdx,1);

uniqueClusters = cellfun(@(c) unique(c), num2cell(clustIdx, 2),'un',0);
uniqueClusters = cellfun(@(c) c(~isnan(c)),uniqueClusters,'un',0);
nClusts = cellfun(@length,uniqueClusters);

clust_change_idx = find(diff(nClusts) ~= 0);

change_k = 1;
for t_k = clust_change_idx(1:end-1)'
    
    current_clust_idxs = clustIdx(t_k:t_k+1,:);
    
    current_clust_nums = unique(current_clust_idxs(1,:));
    next_clust_nums = unique(current_clust_idxs(2,:));
    
    current_clust_nums = current_clust_nums(~isnan(current_clust_nums));
    next_clust_nums = next_clust_nums(~isnan(next_clust_nums));
    
    dead_clusters = [dead_clusters current_clust_nums(~ismember(current_clust_nums,next_clust_nums))];
    if any(ismember(next_clust_nums,dead_clusters))
        clusters2replace = next_clust_nums(ismember(next_clust_nums,dead_clusters));
        for clustNum = clusters2replace
            
            idx2replace = t_k+1:nT;
            clustIdx_new = clustIdx(idx2replace,:);
            clustIdx_new(clustIdx_new == clustNum) = new_clust_num;
            clustIdx(idx2replace,:) = clustIdx_new;
            
%             cluster_duration_idx = find(cellfun(@(c) ~ismember(clustNum,c),uniqueClusters(clust_change_idx(change_k:end)+1)),1,'first');
%             cluster_duration = clust_change_idx(cluster_duration_idx + change_k - 1);
%             cluster_duration = t_k+1:cluster_duration;
%             
%             clustIdx_new = clustIdx(cluster_duration,:);
%             clustIdx_new(clustIdx_new == clustNum) = new_clust_num;
%             
%             clustIdx(cluster_duration,:) = clustIdx_new;
            
            new_clust_num = new_clust_num + 1;
        end
    end
    
    change_k = change_k + 1;
    
end