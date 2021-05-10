function clustIdx = find_co_clustered_bats(pd,expDate)

batPos = pd.batPos(expDate);
batPos = batPos.values;
batPos = cat(3,batPos{:});
batPos = permute(batPos,[3 2 1]);

clustThresh = 15;
nT = size(batPos,3);
nBat = size(batPos,1);
clustIdx = nan(nT,nBat);
nClust = zeros(nT,1);

for t_k = 1:nT
    idx = ~any(isnan(batPos(:,:,t_k)),2);
    if any(idx)
        z = linkage(batPos(idx,:,t_k));
        clustIdx(t_k,idx) = cluster(z,'cutoff',clustThresh,'criterion','distance');
        nClust(t_k) = length(unique(clustIdx(t_k,:)));
    end
end

clustIdx = reassign_clusters(clustIdx,batPos);

end
