function [tp, exp] = dropExtra(tp, exp, cutoff, nkeep)

% this function drops the extra datapoints that are collected for a gene 
% after final expression is reached.
%
% tp        time points
% exp       expression data for one gene at each time point
% cutoff    given as a fraction of the range of exp, used to decide 
%           whether/when final expression is reached. Use higher values 
%           for noisier data. Default 0.05.
% nkeep     how many data points of final expression is kept. Default 2.
%
% Chen Chen. Last update: 2024-08-04
% Rosemary Yu. Last update: 2024-08-08

if nargin < 3
    cutoff = 0.05;
end
if nargin < 4
    nkeep = 2;
end

range = max (exp) - min (exp);
rangep = range * cutoff;
minrange = exp(end) - rangep;
maxrange = exp(end) + rangep;
for i = 1: length (exp)
    if minrange < exp(i) && maxrange > exp(i)
        sign (:,i) = 0;
    else
        sign (:,i) = 1;
    end
end
indexe = width (sign) - 1;
tv = sign (:, indexe);
while tv == 0
    indexe = indexe - 1;
    tv = sign (:, indexe);
end
tpo = tp;
if indexe <= (length(exp) - nkeep -1)
    lsign = indexe + nkeep;
    %extrax = tpo (:,[lsign end]);
    %extray = plotdata2 (g, [lsign end]);
    tp = tpo (:, 1:lsign);
    exp = exp(:, 1:lsign);
end
    
end