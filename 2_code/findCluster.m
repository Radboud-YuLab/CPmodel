function [cluster] = findCluster(desired_cluster, out_opt, cluster_list)
% find genes in a given pre-assigned cluster and pull out all information in out_opt.
% 
% cluster           output data for selected cluster
% desired_cluster   number or name of cluster you want to pull out, can be multiple
%                   in the format of [X; Y]
% out_opt           an m x 6 cell array collecting the second output of the
%                   getOptCP function for m number of genes. Should contain
%                   the following:
%                   col 1: gene name
%                   col 2: number of CPs
%                   col 3: array of CPs. Each row is one CP in (x,y)
%                   coordinates.
%                   col 4: RSS
%                   col 5: nCP-penalized RSS
%                   col 6: R squared
% cluster_list      the list of clusters assigned to each gene according to
%                   GO term analysis by Pommerenke 2012. 
%
% Chen Chen. Last update: 2024-09-11

idx = [];
for i = 1:height (desired_cluster)
    idx_temp = find (cluster_list == desired_cluster (i,:));
    idx = cat (1, idx, idx_temp);
end
cluster = out_opt (idx,:);
end