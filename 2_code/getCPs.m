function output = getCPs(BS, maxnCP, linCutoff, lambda, graphPath, gene, g, gdata)

% main function to fit CP models
%
% BS            the (bootstrapped) data to fit. Given as an n x 2 array: 
%               first col is time point, second col is gene expression.
% maxnCP        max number of CPs to fit. Must be greater than 2 and less 
%               than the number of data points.
% linCutoff     the linear correlation cutoff used to decide whether to
%               fit >2 CPs.
% lambda        the xy-scaling factor used when estimating tau. 
% graphPath     Optional graph parameter. The path to save graphical 
%               outputs. If not given, then graphs are not saved.
% gene          Optional graph parameter. The gene name, as a string.
% g             Optional graph parameter. The number used to order genes.
% gdata         Optional graph parameter. The experimental (non-bootstrapped) 
%               data. If graphs are saved and gdata is given, then gdata is 
%               plotted in a different color.
% 
% output        col 1: number of CPs
%               col 2: array of CPs. Each row is one CP in (x,y)
%               coordinates.
%
% Chen Chen. Last update: 2024-08-04
% Rosemary Yu. Last update: 2024-08-12

if maxnCP < 2
    maxnCP = 2;
end
if maxnCP > length(BS)
    maxnCP = length(BS);
end


% main
output = {};

% for nCP=2, just pull out first and last datapoints
CPs (1,:) = BS (1,:);
CPs (2,:) = BS (end,:);

output (1,1) = num2cell (2);
output (1,2) = {CPs};
%output (1,3) = num2cell (rss);

% check linCutoff to decide whether to fit nCP>2
if maxnCP > 2
    Rtest = corrcoef (BS (:,1), BS (:,2));
    Rtest = abs (Rtest (1,2));
end

% fit nCP>2
if maxnCP > 2 && Rtest < linCutoff
   for c = 2:(maxnCP-1) %c = nCP - 1
       scale_intv = [0, max(BS (:,2)) - min(BS (:,2))] .* lambda; 
       tau = getTau(BS, true, scale_intv);
       B = bernstein_basis (c, tau);
       dataxval = BS (:,1) - (B(:,1).*BS (1,1) + B(:,end).*BS(end,1));
       datayval = BS (:,2) - (B(:,1).*BS(1,2)+ B (:,end).*BS(end,2));
       B (:,1) = [];
       B (:,end) = [];
       
       [bx, ~, ~] = regress (dataxval, B);
       [by, ~, ry] = regress (datayval, B);
       CPs = cat (2, bx, by);
       CPs = cat (1, BS (1,:), CPs, BS (end,:));
       output (c,1) = num2cell (c+1);
       %rss = sum (ry.^2); 
       output (c,2) = {CPs};
       %output (c,3) = num2cell (rss);

   end
end

% if graphPath is given, draw & save graphics
if exist('graphPath', 'var') 
    %check variables
    if ~exist('gene', 'var')
        gene = 'gene_x';        
    end
    gene = convertStringsToChars(gene);
    if ~exist('g', 'var')
        g = 0;        
    end
    
    %make directory
    currentPath = pwd;
    if ~isfolder(graphPath)
        mkdir(graphPath)
    end
    cd (graphPath)
    gfolder = [num2str(g), '_', gene];
    mkdir(gfolder)
    cd (gfolder)
    
    %draw & save graphs
    for ngraphs = 1:height(output)
        P = graph_points(output{ngraphs, 1}, output{ngraphs, 2});
        gcf = figure('visible','off');
        plot(P(:,1), P(:,2));
        hold on
        if exist('gdata', 'var') 
            scatter(BS(:,1), BS(:,2), 10, "black", "o");
            scatter(gdata(:,1), gdata(:,2), "blue", "filled", "o");
        else
            scatter(BS(:,1), BS(:,2), 10, "blue", "o");
        end
        gtitle = [gene, '_', num2str(output{ngraphs, 1}), '_CPs'];
        title(gtitle, 'Interpreter', 'none')
        xlabel ('time')
        ylabel ('expression')
        hold off
        saveas(gcf,gtitle,'png')
    end
    cd (currentPath)
end
end

% local
function B = bernstein_basis(n, tau)
    % Compute Bernstein basis polynomials up to degree n at point tau
    % n = number of control points - 1
    
    % Initialize matrix to store the basis polynomials
    B = zeros(length (tau), n+1);
    
    % Calculate the basis polynomials
    for i = 0:n
        B(:, i+1) = nchoosek(n, i) * tau.^i .* (1 - tau).^(n - i);
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