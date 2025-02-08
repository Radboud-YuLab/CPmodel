function [qu_genes, expr_wt, expr_mut] = readQu

% read dataset from Qu 2018.
% paper DOI: 10.1016/j.celrep.2018.11.039 
% 
% qu_genes      gene names in the Qu 2018 dataset 
% expr_wt       wildtype data (time course days 0-7)
% expr_mut      mutant data (3 different mutants, in the order 
%               [R204W R279H R304W])
% 
% Chen Chen. Last update: 2025-02-06
% Rosemary Yu. Last update: 2025-02-07


% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 12);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["genes", "wt_d0", "wt_d1", "wt_d2", "wt_d3", "wt_d4", "wt_d5", "wt_d6", "wt_d7","R204W_d7", "R279H_d7", "R304W_d7"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "genes", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "genes", "EmptyFieldRule", "auto");

% Import the data
qu_data = readtable("1_raw_data\Qu_2018_wt_mut_clean.txt", opts);

% pull out gene names & data
genes_all = table2array(qu_data(:,1));
expr_wt_all =  table2array(qu_data(:,2:9));
expr_mut = table2array(qu_data(:,10:12));


% Filter for only DE genes
DE_wt = sum (expr_wt_all >= 1 | expr_wt_all <= -1, 2);
idxDE = find (DE_wt~=0);
expr_de = expr_wt_all (idxDE,:);
qu_genes = genes_all (idxDE,:);
expr_mut = expr_mut (idxDE,:);

% Filter for genes that are DE at at least 2 timepoints 
filter_logFC = sum(expr_de~=0,2);
idx_filter_logFC = find(filter_logFC >= 2);
expr_wt = expr_de (idx_filter_logFC,:);
qu_genes = qu_genes (idx_filter_logFC,:);
expr_mut = expr_mut (idx_filter_logFC,:);


end