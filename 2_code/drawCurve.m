function [gcf] = drawCurve(out_opt, pom_days, pom_log2_fc, gene, inputmode, path)
% produces a figure with the data plotted, as well as the optimal CPs and
% the resulting curve if applicable. 
%   
% gcf               figure output
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
% pom_days          time points
% pom_log2_fc       expression data for the genes at each time point
% gene              gene name or number to be plotted
% inputmode         1 for gene name
%                   2 for gene number
% path              the path to save the graphical outputs
%
% Chen Chen. Last update: 2024-09-06

for i = 1:height (gene)
    if inputmode == 1
        idx = regexp (gene (i,:), out_opt (:,1));
        idx = find (~cellfun (@isempty, idx));
        geneplot = out_opt (idx,:);
        log2_fc = pom_log2_fc (idx,:);
        name = gene
    elseif inputmode == 2
        geneplot = out_opt (gene (i,:),:);
        log2_fc = pom_log2_fc (gene (i,:),:);
        name = out_opt {gene (i,:),1};
    end

    data = transpose (cat (1, pom_days, log2_fc));
    gcf = figure('visible','off');
    hold on
    if ~isempty (geneplot {1,2})
        P = graph_points (geneplot {:,2}, geneplot {:,3});
        plot(P(:,1), P(:,2),"Color",[0 0.64 0.93],"LineWidth",3.0);
        scatter (geneplot {:,3} (:,1), geneplot {:,3} (:,2), 40, [0 0.64 0.93],"filled");
        ncp = num2str (geneplot {:,2});
        gtitle = [name, ' ',  ncp, ' CPs'];
        gtitle2 = [name, '_',  ncp, '_CPs'];
    else
        gtitle = [name];
        gtitle2 = gtitle;
    end
    scatter (data (:,1), data (:,2), 10,"black","filled");
    title (gtitle)
    xlabel ('time')
    ylabel ('expression')
    if exist('path', 'var')
        %make directory
        currentPath = pwd;
        if ~isfolder(path)
            mkdir(path)
        end
        cd (path)
        saveas(gcf,gtitle2,'png')
        cd (currentPath)
    end
end
end


function P = graph_points(n, cp)
    % generate points from CPs for graphing
    % n = number of CPs
    
    n1 = n-1;

    % calculate the binomial coefficient (x!/(y!(x-y)!)) 
    for i = 0:1:n1
        binomial_coefficient(i+1) = factorial(n1)/(factorial(i)*factorial(n1-i));   
    end

    basis = [];
    UB = [];
    for tau = 0:0.01:1
        for k = 1:n
            UB(k) = binomial_coefficient(k)*((1-tau)^(n-k))*(tau^(k-1));
        end
        basis = cat(1,basis,UB);    
    end

    P = basis*cp;

end