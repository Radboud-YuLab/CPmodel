function BS = bootTimeSeries(tp, exp, nsample, tau)

% Bootstrap a time-series dataset. Closer timepoints have higher weights.
%
% tp        time points
% exp       expression data for one gene at each time point
% nsample   how many points to sample. Default 100.
% tau       if given, return data is linearly spaced on tau. If not
%           given, return data is linearly spaced on tp.
%
% Chen Chen. Last update: 2024-08-04
% Rosemary Yu. Last update: 2024-08-09

if nargin < 3
    nsample = 100;
end

if exist('tau', 'var')
    %if taus is given, linspace on tau
    tau_ls = linspace(1e-5,1-1e-5,nsample).';
    BS = interp1(tau,tp,tau_ls);
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
else
    %if tau is not given, linspace on time
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