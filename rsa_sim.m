
function [rwith, rbetw] = rsa_sim(pattern,corrmat)

% get within / between condition similarities

[rwith, rbetw] = deal([]);

%get pairs for between similiarity in both runs

pattern.pairs = combnk(unique(pattern.tLabel),2);

%get within-similiarity

for s = 1:numel(unique(pattern.tLabel))
    
    clear idx; idx = find(pattern.tLabel == s);
    rwith = [rwith, nanmean(reshape(corrmat.fisher(min(idx):max(idx),min(idx):max(idx)),1,numel(idx)^2),2)];
    
end

%get between-similiarity

for s = 1:size(pattern.pairs,1)
    
    clear idx; idx = find(pattern.tLabel == pattern.pairs(s,1));
    clear iidx; iidx = find(pattern.tLabel == pattern.pairs(s,2));
    if min(idx) > min(iidx); clear idx; idx = iidx; clear iidx; iidx = find(pattern.tLabel == pattern.pairs(s,1)); end
    rbetw = [rbetw, nanmean(reshape(corrmat.fisher(min(idx):max(idx),min(iidx):max(iidx)),1,numel(idx)*numel(iidx)),2)];
    
end
