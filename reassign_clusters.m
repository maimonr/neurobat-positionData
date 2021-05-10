function clustIdx = reassign_clusters(clustIdx,batPos)

clustIdx_orig = clustIdx;

nT = size(clustIdx,1);
nBat = size(clustIdx,2);

clust_change_idx = false(1,nT);
for t_k = 1:nT-1
    clust_change_idx(t_k) = ~isequaln(clustIdx_orig(t_k,:),clustIdx_orig(t_k+1,:));
end

clust_change_idx = [find(clust_change_idx) nT];

change_k = 1;
for t_k = clust_change_idx(1:end-1)
    current_clust_nums = clustIdx(t_k:t_k+1,:);
    
    if all(isnan(current_clust_nums(1,:)))
        change_k = change_k + 1;
        continue
    end
    
    cluster_nums = unique(current_clust_nums(:));
    cluster_nums = cluster_nums(~isnan(cluster_nums))';
    nClust = length(cluster_nums);
    
    p = cell(1,2);
    for step_k = 1:2
        if step_k == 1
            p{step_k} = inf(nClust,2);
        else
            p{step_k} = nan(nClust,2);
        end
        clust_k = 1;
        for cluster_num = cluster_nums
            if any(current_clust_nums(step_k,:) == cluster_num)
                p{step_k}(clust_k,:) = nanmean(batPos(current_clust_nums(step_k,:) == cluster_num,:,t_k + step_k - 1),1);
            end
            clust_k = clust_k + 1;
        end
    end
    
    clustIdx_new = nan(1,nBat);
    
    clustDist = squareform(pdist(vertcat(p{:})));
    clustDist = clustDist(1:nClust,nClust+1:end);
    [~,sortIdx] = sort(min(clustDist));
    
    clust_k = 1;
    while any(isnan(clustIdx_new)) && clust_k <= length(cluster_nums) && sortIdx(clust_k) <= size(p{2},1)
        current_cluster_num_idx = sortIdx(clust_k);
        cluster2reassign = cluster_nums(current_cluster_num_idx);
        d = vecnorm(p{2}(current_cluster_num_idx,:) - p{1},2,2);
        if all(isnan(d))
            break
        end
        [~,minIdx] = min(d);
        clustIdx_new(current_clust_nums(2,:) == cluster2reassign) = cluster_nums(minIdx);
        p{1}(minIdx,:) = NaN;
        clust_k = clust_k + 1;
    end
    
    assert(length(unique(current_clust_nums(2,:))) == length(unique(clustIdx_new)))
    
    rows2change = t_k+1:clust_change_idx(change_k+1);
    clustIdx(rows2change,:) = repmat(clustIdx_new,length(rows2change),1);
    change_k = change_k + 1;
end

end