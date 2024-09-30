function pred_exp = predUhlitzCYHX(getOptCPoutput, predTime, timeCutoff)

% calculates the predicted effect (on gene expression) of perturbing a 
% completion process, based on the fitted CPs of the unperturbed process.
%
% getOptCPoutput   an m x 3 cell array of the second output of the
%                  getOptCP function for m number of genes. Should contain
%                  the following:
%                  col 1: gene name
%                  col 2: number of CPs
%                  col 3: array of CPs. Each row is one CP in (x,y) coordinates.
%                  additional columns are allowed but will not be used.
% predTime         an array
% timeCutoff       the estimated time for the state transition that is
%                  blocked by CHYX treatment
% pred_exp         an m x n array of the predicted expression levels. m is the 
%                  number of genes (from out_opt) and n is the number of time 
%                  points (from predTime).
%
% Rosemary Yu. Last update: 2024-08-22

nGenes = height(getOptCPoutput);
pred_exp = NaN(nGenes,length(predTime));

for g = 1:nGenes
    gCP = cell2mat(getOptCPoutput(g,3));
    % if there are no optCPs, skip
    if isempty(gCP)
        continue
    end
    
    % if there are intCPs < timeCutoff, pull out all CPs < timeCutoff
    if gCP(2,1) < timeCutoff
        gCP = gCP( gCP(:,1) < timeCutoff, : );
    end
        
    for idx = 1:length(predTime)
        tp = predTime(idx);
        %fill in predTime >= last CP 
        if tp >= gCP(end,1)
            pred_exp(g,idx) = gCP(end,2);
        else
            %for predTime < last CP, fit from CPs (see getOptCP)
            %calculate r from tp
            b = gCP(:,1);
            poly = changeBasis(b);
            poly(end) = poly(end) - tp;
            r = roots(poly);
            r = real(r(imag(r) == 0));
            r = r(round(r,4)>=0 & round(r,4) <=1);

            % use r to calculate predicted expression
            pred_exp(g,idx) = calculate_y(length(gCP(:,2)), gCP(:,2), r);
        end
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