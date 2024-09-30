function [outAll, outOpt] = getOptCP(gene, gdata, getCPoutput, optType, optCutoff) 

% Select the optimal number of CPs for a gene time-series data. 
%
% gene          gene name
% gdata         an n x 2 array of the time-series data (experimental, 
%               not boot-strapped). First column is the x-variable  
%               (time), second column is the y-variable (gene 
%               expression).
% getCPoutput   output of the getCPs function. Should contain the
%               following:
%               col 1: number of CPs
%               col 2: array of CPs. Each row is one CP in (x,y)
%               coordinates.
% optType       what is used to select for optimal number of CPs.
%               Options are: 
%               1 - min penalized RSS (default)
%               2 - cutoff of scaled RSS
%               3 - cutoff of R2
% OptCutoff     If optType is 2 or 3, then this is the cutoff to be used. 
% 
% outAll        an array with goodness-of-fit measurements for all fitted 
%               models.
%               col 1: gene name
%               col 2: number of CPs
%               col 3: array of CPs. Each row is one CP in (x,y)
%               coordinates.
%               col 4: RSS
%               col 5: nCP-penalized RSS
%               col 6: R squared
% outOpt        Contains the same columns as outAll but only for the model
%               with the optimal number of CPs. 
%
% Chen Chen. Last update: 2024-08-04
% Rosemary Yu. Last update: 2024-08-13

if nargin < 4
    optType = 1;
end
if optType ~= 1 && ~exist('optCutoff', 'var')
    error('cutoff value needed')
end

outAll = {};

%calculate RSS, pRSS, and R-squared
for j = 1:height(getCPoutput)
    nCP = cell2mat(getCPoutput(j,1));
    b = getCPoutput{j,2}(:,1);
    if ~all(diff(b)>0)
        % if CPs on x (time) are not monotonously increasing, model is dropped
        continue 
    end
    c = changeBasis(b);
    y_hat = zeros(height(gdata),1);

    for k = 1:height(gdata)
        % calculate r from x-value
        poly = c;
        poly(end) = poly(end) - gdata(k,1);
        r = roots(poly);
        r = real(r(imag(r) == 0));
        r = r(round(r,4)>=0 & round(r,4) <=1);

        % use r to calculate fitted y-value (y_hat)
        P = calculate_y(nCP, getCPoutput{j,2}(:,2), r);
        y_hat(k) = P;
    end
    
    %calculate RSS and nCP-penalized RSS
    rss = sum((gdata(:,2) - y_hat).^2);
    sp_rss = rss * (nCP-1)^2;
    
    %calculate R squared
    sst = sum((gdata(:,2) - mean(gdata(:,2))).^2);
    rsq = 1-(rss/sst);
    if rsq < 0
        rsq = 0;
    end
    
    outAll(j,1) = cellstr(gene);
    outAll(j,2) = num2cell(nCP);
    outAll(j,3) = {getCPoutput{j,2}};
    outAll(j,4) = num2cell(rss);
    outAll(j,5) = num2cell(sp_rss);
    outAll(j,6) = num2cell(rsq);
    
end

outAll(all(cellfun(@isempty, outAll),2),:) = [];

if optType == 1 
    % find nCP with min penalized RSS
    minval = min(cell2mat(outAll(:,5)));
    outOpt = outAll(cell2mat(outAll(:,5))==minval,:); 
end

if optType == 2
    %find nCP that meets cutoff of scaled RSS
    bfd = normalize(cell2mat(outAll(:,4)), "norm", Inf);
    nocp = find (bfd < optCutoff); 
    
    %if nocp is empty, check if there are fitted CP>2 models
    if isempty (nocp) 
        if height(getCPoutput)==1
            outOpt = outAll(1,:);
        else
            outOpt = outAll(1,:);
            outOpt(1,2:end) = {[],[],[],[],[]};
        end
    %if nocp exists, return the first nCP that meets the cutoff
    else
        outOpt = outAll(nocp(1,1),:);
    end
end

if optType ==3
    %find nCP that meets cutoff of R2
    R2 = cell2mat(outAll(:,6));
    nocp = find(R2 > optCutoff);
    
    %if nocp is empty, keep only the gene name & return nothing
    if isempty(nocp)
        outOpt = outAll(1,:);
        outOpt(1,2:end) = {[],[],[],[],[]};
    %if nocp exists, return the first nCP that meets the cutoff
    else
        outOpt = outAll(nocp(1,1),:);
    end
end

end

% local
function c = changeBasis(b)
% convert Bernstein basis control points to power basis coefficients.

c = zeros(length(b),1);
c(1) = b(1);
res = b;
for i = 2:length(b)
    u = res(2:end);
    v = res(1:(end-1));
    res = u - v;
    c(i) = res(1) * nchoosek(length(b)-1, i-1);
end
c = flip(c);

end

function P = calculate_y(n, cp, tau)
    % calculate a point by tau and CPs. n = number of CPs
    
    n1 = n-1;

    % calculate the binomial coefficient (x!/(y!(x-y)!)) 
    for i = 0:1:n1
        binomial_coefficient(i+1) = factorial(n1)/(factorial(i)*factorial(n1-i));   
    end

    basis = [];
    UB = [];
    for k = 1:n
        UB(k) = binomial_coefficient(k)*((1-tau)^(n-k))*(tau^(k-1));
    end
    basis = cat(1,basis,UB);    

    P = basis*cp;

end