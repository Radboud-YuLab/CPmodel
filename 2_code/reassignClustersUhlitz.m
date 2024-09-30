function [IEG, DEG, SRG, cluster_list] = reassignClustersUhlitz
%Clusters the genes based on an arbitrary cutoff value
%   IEG = immediate early genes
%   DEG = delayed early genes
%   SRG = secondary response genes
%
% Rosemary Yu. Last update: 2024-09-11
% Chen Chen. Last update: 2024-09-14

% Determine SRGs: genes that do not increase expression in the CYHX dataset
[CYHX_names, CYHX_log2fc, ~] = readUhlitzCYHX;
for i = 1:height (CYHX_log2fc)
    expr = find ((CYHX_log2fc (i,:) <= 1) == 0); %arbitrary cutoff at log2fc of 1
    if ~isempty (expr)
        idxSRG (i,:) = 0;
    else
        idxSRG (i,:) = 1;
    end
end
idxSRG = find (idxSRG == 1);
CYHX_names_SRG = CYHX_names (idxSRG,:);
CYHX_log2fc = CYHX_log2fc (idxSRG,:);

for i = 1:height (idxSRG)
    cluster_list (idxSRG (i,:), :) = "SRG";
end

SRG = table (CYHX_names_SRG, CYHX_log2fc (:,1), CYHX_log2fc (:,2), CYHX_log2fc (:,3));
SRG.Properties.VariableNames = ["gene name", "T1h_CYHX_1h", "T1h_CYHX_2h", "T1h_CYHX_4h"];

% Process data from 4OHT experiment
[OHT_names, ~, OHT_log2fc] = readUhlitzTimeSeries;
OHT_log2fc (idxSRG,:) = [];
OHT_names (idxSRG, :) = [];

% Determine IEGs: genes that increase expression within 1 h
IEGtest = OHT_log2fc;

IEGtest = IEGtest (:,2:3); 
for i = 1:height (IEGtest)
    expr = find ((IEGtest (i,:) < 1.0) == 0); %arbitrary cutoff
    if ~isempty (expr)
        idxIEG (i,:) = 1;
    else
        idxIEG (i,:) = 0;
    end
end
idxIEG = find (idxIEG == 1);
IEG_names = OHT_names (idxIEG,:);
IEG_log2fc = OHT_log2fc (idxIEG,:);

for i = 1:height (IEG_names)
    trueidxIEG = find (strcmp (IEG_names (i,:), CYHX_names));
    cluster_list (trueidxIEG, :) = "IEG";
end

IEG = table (IEG_names, IEG_log2fc (:,1), IEG_log2fc (:,2), IEG_log2fc (:,3), IEG_log2fc (:,4), IEG_log2fc (:,5), IEG_log2fc (:,6), IEG_log2fc (:,7), IEG_log2fc (:,8));
IEG.Properties.VariableNames = ["gene name", "T05h", "T1h", "T2h", "T3h", "T4h", "T6h", "T8h", "T10h"];

% Determine DEG: all genes that are not considered IEG or SRG
DEG_names = OHT_names;
DEG_log2fc = OHT_log2fc;
DEG_log2fc (idxIEG,:) = [];
DEG_names (idxIEG, :) = [];

for i = 1:height (DEG_names)
    trueidxDEG = find (strcmp (DEG_names (i,:), CYHX_names));
    cluster_list (trueidxDEG, :) = "DEG";
end

DEG = table (DEG_names, DEG_log2fc (:,1), DEG_log2fc (:,2), DEG_log2fc (:,3), DEG_log2fc (:,4), DEG_log2fc (:,5), DEG_log2fc (:,6), DEG_log2fc (:,7), DEG_log2fc (:,8));
DEG.Properties.VariableNames = ["gene name", "T05h", "T1h", "T2h", "T3h", "T4h", "T6h", "T8h", "T10h"];

end