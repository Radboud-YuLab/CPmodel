function [CYHX_names, CYHX_log2fc, CYHX_clusters] = readUhlitzCYHX

% this function reads and processes the Uhlitz CYHX treatment dataset. 
% The dataset is accessible at: 
% https://doi.org/10.1371/journal.pone.0041169.s003
% The paper is accessble at:
% https://doi.org/10.15252/msb.20177554

% Chen Chen. Last update: 2024-09-11
% Rosemary Yu. Last update: 2024-09-12


% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 8);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["ID_REF", "SYMBOL", "cluster", "Var4", "Var5", "T1h_CYHX_1h", "T2h_CYHX_2h", "T4h_CYHX_4h"];
opts.SelectedVariableNames = ["ID_REF", "SYMBOL", "cluster", "T1h_CYHX_1h", "T2h_CYHX_2h", "T4h_CYHX_4h"];
opts.VariableTypes = ["double", "string", "string", "string", "string", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["SYMBOL", "Var4", "Var5"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["SYMBOL", "cluster", "Var4", "Var5"], "EmptyFieldRule", "auto");

% Import the data
Uhlitzlog2FCCYHX = readtable("1_raw_data\Uhlitz_log2FC_CYHX.csv", opts);

% Clear temporary variables
clear opts

CYHX_log2fc = table2array (Uhlitzlog2FCCYHX (2:end, 4:6));
CYHX_names = table2array (Uhlitzlog2FCCYHX (2:end,2));
CYHX_clusters = table2array (Uhlitzlog2FCCYHX (2:end, 3));
end
