function clustIdx = clean_clusters(clustIdx)

nBat = size(clustIdx,2);

video_fs = 20;
min_bat_in_cluster = 3;
min_cluster_duration = 10;
min_cluster_frames = video_fs*min_cluster_duration;

all_clust_nums = unique(clustIdx(~isnan(clustIdx)));

for cNum = all_clust_nums'
    n_bats_in_clust = sum(clustIdx == cNum,2);
    discardIdx = repmat(n_bats_in_clust < min_bat_in_cluster,1,nBat) & clustIdx == cNum;
    clustIdx(discardIdx) = 0;
end

% n_bats_in_clust = cellfun(@(cNum) max(sum(clustIdx == cNum,2)),num2cell(all_clust_nums));
% clusters_to_discard = all_clust_nums(n_bats_in_clust < min_bat_in_cluster);
% discardIdx = ismember(clustIdx,clusters_to_discard);
% clustIdx(discardIdx) = 0;

all_clust_nums = unique(clustIdx(~isnan(clustIdx)));
all_clust_nums = all_clust_nums(all_clust_nums~=0);

clusterDuration = cellfun(@(cNum) sum(any(clustIdx == cNum,2)),num2cell(all_clust_nums));
clusters_to_discard = all_clust_nums(clusterDuration < min_cluster_frames);
discardIdx = ismember(clustIdx,clusters_to_discard);
clustIdx(discardIdx) = 0;

end