function [pts,idx] = remove_duplicate_points(pts)
[~,I,~] = unique(round(pts), 'rows', 'first');
ixDupRows = setdiff(1:size(pts,1), I);
idx = 1:size(pts,1);
pts(ixDupRows,:) = [];
idx(ixDupRows) = [];