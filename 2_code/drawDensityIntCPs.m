function drawDensityIntCPs(out_opt, timeRange, intOptions, bandwidth, graphPath, graphName)

% draw density plot of the intermediate CPs in the time range of the
% timeseries dataset.
%
% out_opt          an m x 6 cell array collecting the second output of the
%                  getOptCP function for m number of genes. Should contain
%                  the following:
%                  col 1: gene name
%                  col 2: number of CPs
%                  col 3: array of CPs. Each row is one CP in (x,y)
%                  coordinates.
%                  col 4: RSS
%                  col 5: nCP-penalized RSS
%                  col 6: R squared
% timeRange        an array of time range of the dataset. Can be just the
%                  start and end time (e.g. [0, 60]) or all timepoints of
%                  the timeseries (e.g. [0,1,2,3,...,60])
% intOptions       What is plotted in the density plot. Options are:
%                  1 - first positive intermediate CPs (default)
%                  2 - all intermediate CPs
% bandwidth        bandwith used for kdensity. Default 0.1.
% graphPath        Optional graph parameter. The path to save graphical
%                  outputs. If not given, then graphs are not saved.
% graphName        Optional to assign the gene name to the graph. If not
%                  given, the gene name will not be included.
%
% Rosemary Yu. Last update: 2024-08-15
% Chen Chen. Last update: 2024-09-14


if nargin < 4
    bandwidth = 0.1;
end
if nargin < 3
    intOptions = 1;
end

plotdata = [];

for i = 1:height(out_opt)
    nIS = cell2mat(out_opt(i,2))-2;
    if ~isempty(nIS) && nIS > 0
        geneIS = cell2mat(out_opt(i,3));
        if intOptions == 1
            %first IS of each gene
            idxnz = find (geneIS (:,2) > 0);
            if ~isempty (idxnz) && geneIS (idxnz (1,1),1) < timeRange (end)
                geneIS = geneIS(idxnz (1,1), 1);
                plotdata = [plotdata; geneIS];
                if ~exist('graphName', 'var')
                    gtitle = 'density_first_intermediate_CPs';
                else
                    gtitle = [graphName, '_density_first_intermediate_CPs'];
                end
            end
        elseif intOptions == 2
            %all IS's of each gene
            geneIS = geneIS(2:(end-1),1);
            plotdata = [plotdata; geneIS];
            if ~exist('graphName', 'var')
                gtitle = 'density_all_intermediate_CPs';
            else
                gtitle = [graphName, '_density_all_intermediate_CPs'];
            end
        end

    end
end

%make plot
[f,xi] = ksdensity(plotdata, 'Support', [timeRange(1) timeRange(end)], 'Bandwidth', bandwidth);
gcf = figure;
plot(xi,f)
xlim([timeRange(1) timeRange(end)])
%title(gtitle, 'Interpreter', 'none')
xlabel ('time');
ylabel ('kernel density of CPs');


% if graphPath is given, save graphics
if exist('graphPath', 'var')
    %make directory
    currentPath = pwd;
    if ~isfolder(graphPath)
        mkdir(graphPath)
    end
    cd (graphPath)

    saveas(gcf,gtitle,'png')
    cd (currentPath)
end


end