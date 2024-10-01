% main script for data analysis and visualization in the manuscript
% "A multi-step completion process model of cell plasticity".
% The folders 1_raw_data and 2_code must already exist with the appropriate files.
% Other folders and files will be written (or over-written) by the script.

% Chen Chen. Last update: 2024-09-16
% Rosemary Yu. Last update: 2024-09-25

%% Pommerenke 2012 dataset, viral infection-induced immune response
close all
clear

% read in the Pommerenke 2012 dataset
pom_days = [0, 1, 2, 3, 5, 8, 10, 14, 18, 22, 26, 30, 40, 60];
[pom_gene_names, pom_log2_fc, pom_clusters] = readPommerenke;

% tune/set data-specific parameters 
data = pom_log2_fc;
time_points = pom_days;
gene_names = pom_gene_names;
linCutoff = 0.95;
lambda = 4;
R2Cutoff = 0.6;

% set up output structures
out_all = struct();
out_opt = {};
out_opt_default = {};  

% fit CP model for each gene
tic

for g = 1:length(data)
    gene = gene_names(g);
    tp = time_points;
    exp = data(g, :);
    
    % drop "extra" datapoints after final expression reached
    [tp, exp] = dropExtra(tp, exp, 0.1, 2); 
    
    % bootstrap
    gdata = horzcat(tp.', exp.');
    interval_x = [0, max(exp) - min(exp)]; 
    tau = getTau(gdata, true, interval_x);
    BS = bootTimeSeries(tp, exp, 100, tau); 
    
    % get CPs
    maxnCP = floor((length(tp)-2)/3)+2;
    output = getCPs(BS, maxnCP, linCutoff, lambda);
    
    % optional: to save graphics (note that run time will be much longer)
    %graphPath = '4_processed_data\pom_getCPs_graphOutput';
    %output = getCPs(BS, maxnCP, linCutoff, lambda, graphPath, gene, g, gdata);
    
    % get optimal number of CPs using R2Cutoff. Set to relatively low (0.6)
    % since data is noisy
    [gAll, gOpt] = getOptCP(gene, gdata, output, 3, R2Cutoff);
    out_opt(g,1:6) = gOpt;
    
    % alternatively can use getOptCP with default parameters: 
    [~, out_opt_default(g,1:6)] = getOptCP(gene, gdata, output);
    
    % grab all data
    out_all(g).gene = gene;
    out_all(g).n_CPs = cell2mat(gAll(:,2));
    out_all(g).CP_coordinates = gAll(:,3);
    out_all(g).RSS = cell2mat(gAll(:,4));
    out_all(g).nCP_penalized_RSS = cell2mat(gAll(:,5));
    out_all(g).R2 = cell2mat(gAll(:,6));
    
end

toc

% save outputs as a matlab object
currentPath = pwd;
mkdir('4_processed_data\pom_main')
cd('4_processed_data\pom_main')
save('pom_main.mat', 'out_all', 'out_opt', 'out_opt_default', 'pom_days', 'pom_clusters', 'pom_log2_fc');
cd (currentPath)
% to load the output:
%load('4_processed_data\pom_main\pom_main.mat')

% draw & save kernel density graph
drawDensityIntCPs(out_opt, pom_days, 1, 0.1, '6_results\pom_main')

% kernel density with getOptCP with default parameters:
out_opt_def_R2cutoff = out_opt_default(find(cell2mat(out_opt_default(:,6)) > 0.6),:);
drawDensityIntCPs(out_opt_def_R2cutoff, pom_days, 1, 0.1)

% write output as csv
[nCP,CPloc] = writeOptCP(out_opt, '6_results\pom_main');

% pull out different clusters in Pommerenke 2012 & draw curves
%   Cluster 6 - innate response
%   Cluster 4 - T cell response
%   Cluster 1 and 2 - B cell response
%   Cluster 7 - tissue repair

[cluster6] = findCluster (6, out_opt, pom_clusters);
drawDensityIntCPs(cluster6, pom_days, 1, 0.1, '6_results\pom_main', 'cluster6')
[cluster4] = findCluster (4, out_opt, pom_clusters);
drawDensityIntCPs(cluster4, pom_days, 1, 0.1, '6_results\pom_main', 'cluster4')
[cluster1_2] = findCluster ([1;2], out_opt, pom_clusters);
drawDensityIntCPs(cluster1_2, pom_days, 1, 0.1, '6_results\pom_main', 'cluster1_2')
[cluster7] = findCluster (7, out_opt, pom_clusters);
drawDensityIntCPs(cluster7, pom_days, 1, 0.1, '6_results\pom_main', 'cluster7')

%% Uhlitz 2017 dataset, ERK signaling
close all
clear

% read in Uhlitz 2017 dataset
uhl_hrs = [0,0.5,1,2,3,4,6,8,10];
[uhl_gene_names, ~, uhl_log2_fc] = readUhlitzTimeSeries;

% tune/set data-specific parameters 
data = uhl_log2_fc; 
time_points = uhl_hrs;
gene_names = uhl_gene_names;
linCutoff = 0.99;
lambda = 4.5;
R2Cutoff = 0.85;

% set up output structures
out_all = struct();
out_opt = {};

% fit CP model for each gene
tic

for g = 1:length(data)
    gene = gene_names(g);
    tp = time_points;
    exp = data(g, :);
    
    % drop "extra" datapoints after final expression reached
    [tp, exp] = dropExtra(tp, exp, 0.01, 2); 
    
    % bootstrap
    gdata = horzcat(tp.', exp.');
    interval_x = [0, max(exp) - min(exp)]; 
    tau = getTau(gdata, true, interval_x);
    BS = bootTimeSeries(tp, exp, 100, tau); 
    
    % get CPs
    maxnCP = floor((length(tp)-2)/2)+2;
    output = getCPs(BS, maxnCP, linCutoff, lambda);
    
    % optional: to save graphics (note that run time will be much longer)
    %graphPath = '4_processed_data\uhl_getCPs_graphOutput';
    %output = getCPs(BS, maxnCP, linCutoff, lambda, graphPath, gene, g, gdata); 
    
    % get optimal number of CPs. 
    [gAll, gOpt] = getOptCP(gene, gdata, output, 3, R2Cutoff);
    out_opt(g,1:6) = gOpt;
    
    % grab all data
    out_all(g).gene = gene;
    out_all(g).n_CPs = cell2mat(gAll(:,2));
    out_all(g).CP_coordinates = gAll(:,3);
    out_all(g).RSS = cell2mat(gAll(:,4));
    out_all(g).nCP_penalized_RSS = cell2mat(gAll(:,5));
    out_all(g).R2 = cell2mat(gAll(:,6));
    
end

toc

% Reassign the genes into IEG, DEG and SRG
[IEG, DEG, SRG, uhl_cluster] = reassignClustersUhlitz;

% save outputs as a matlab object
currentPath = pwd;
mkdir('4_processed_data\uhl_main')
cd('4_processed_data\uhl_main')
save('uhl_main.mat', 'out_all', 'out_opt', 'uhl_hrs', 'uhl_cluster');
cd (currentPath)
% to load the object:
%load('4_processed_data\uhl_main\uhl_main.mat')

% draw & save kernel density graph
drawDensityIntCPs(out_opt, uhl_hrs, 1, 0.2, '6_results\uhl_main')

% write output as csv
[nCP,CPloc] = writeOptCP(out_opt, '6_results\uhl_main');

% pull out genes in each cluster & draw curves
[SRG_cluster] = findCluster ('SRG', out_opt, uhl_cluster);
drawDensityIntCPs(SRG_cluster, uhl_hrs, 1, 0.2, '6_results\uhl_main', 'clusterSRG');
[IEG_cluster] = findCluster ('IEG', out_opt, uhl_cluster);
drawDensityIntCPs(IEG_cluster, uhl_hrs, 1, 0.2, '6_results\uhl_main', 'clusterIEG');


%% predict the response of interrupting a state transition
close all
clear

load('4_processed_data\uhl_main\uhl_main.mat')
timeCutoff = 2;
predTime = [1 2 4];
pred_exp = predUhlitzCYHX(out_opt, predTime, timeCutoff);

%calculate R2
[~, CYHX_log2fc, ~] = readUhlitzCYHX;

rss = sum((CYHX_log2fc - pred_exp).^2, 'omitnan');
sst = sum((CYHX_log2fc - mean(CYHX_log2fc)).^2);
rsq = 1 - rss./sst;

%calculate R2 of null model (random permutaion of the CYHX experimental values)
idx = randperm(numel(CYHX_log2fc));
rand_CYHX = CYHX_log2fc(idx(1:numel(pred_exp)));
pred_null = reshape(rand_CYHX, size(pred_exp));

rss_null = sum((CYHX_log2fc - pred_null).^2, 'omitnan');
rsq_null = 1 - rss_null./sst;

%save files
currentPath = pwd;
cd('4_processed_data\uhl_main')
save('uhl_pred_CYHX.mat', 'pred_exp', 'pred_null', 'predTime', 'CYHX_log2fc', 'rsq', 'rsq_null');
cd (currentPath)
% to load the object:
%load('4_processed_data\uhl_main\uhl_pred_CYHX.mat')

T = table(out_opt(:,1), pred_exp(:,1), pred_exp(:,2), pred_exp(:,3), ...
    'VariableNames', {'SYMBOL', 'T1h_CYHX', 'T2h_CYHX', 'T4h_CYHX'});
writetable(T, '6_results\uhl_main\CHYX_predictions.csv')


%% draw figures
close all
clear

mkdir('6_results\figures')
[Figure1d,Figure2a,Figure2b,Figure3] = drawFigures;

disp('all done')
