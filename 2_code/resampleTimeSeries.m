function BS = resampleTimeSeries(tp, exp, nsample)

% Resample a time-series dataset, linearly spaced on the time variable. 
% Closer timepoints have higher weights.
%
% tp        time points
% exp       expression data for one gene at each time point
% nsample   how many points to sample. Default 100.
%
% Chen Chen. Last update: 2024-08-04
% Rosemary Yu. Last update: 2024-08-08

if nargin < 3
    nsample = 100;
end


BS = linspace(0, tp(end), nsample).';
for i = 1:length(BS)
    if ismember (BS(i), tp)
        inx = find (tp == BS(i));
        BS (i, 2) = exp (inx);
    else
        dist = 1./ (abs(BS(i)-tp).^2.5);
        weight = normalize (dist, 'norm', 1);
        bsd = randsample(exp, 1000, true, weight);
        bsd = sum (bsd)/1000;
        BS (i,2) = bsd;
    end
end

end