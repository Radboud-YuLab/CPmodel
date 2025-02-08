function [genes] = findState(statelimit,gene_state_no, out_opt, timeRange, path)
% Pulls out the genes in each intermediate attractor state and puts it into a txt file. 
%
% statelimit:       upper and lower bound of the time value for the 
%                   intermediate attractor state (in [t1 t2] format)
% gene_state_no:    number of the intermediate state
% out_opt:          out_opt as described in getOptCP
% timeRange:        total time range of the data
% path:             pathway in which the txt file should be saved
%
% Chen Chen. Last update: 2025-02-06

%filter out the genes for which it could find an optimal amount of CPs
empty = cellfun (@isempty, out_opt (:,2));
idx = find (empty);
out_opt2 = out_opt;
out_opt2 (idx,:) = [];
list = zeros(height(out_opt2),1);

for i = 1:height (out_opt2)
    geneIS = cell2mat(out_opt2(i,3));
    nIS = cell2mat(out_opt2(i,2))-2;
    if ~isempty(nIS) && nIS > 0
        geneIS = cell2mat(out_opt2(i,3));
        idxnz = find (geneIS (:,2) > 0);
        if ~isempty (idxnz) && geneIS (idxnz (1,1),1) < timeRange (end)
            geneIS = geneIS(idxnz (1,1), 1);
            
            if geneIS > statelimit (1,1) && geneIS < statelimit (1,2);
                list (i,:) = 1;
            else
                list (i,:) = 0;
            end
        end
        %geneIS = geneIS(2,1);
        %if geneIS > statelimit (1,1) && geneIS < statelimit (1,2);
            %list (i,:) = 1;
        %else
            %list (i,:) = 0;
        %end
    else
        list (i,:) = 0;
    end
end
idxlist = find (list == 0);
out_opt2 (idxlist,:) = [];
genes = out_opt2 (:,1);

if exist('path', 'var') 
    %make directory
    currentPath = pwd;
    if ~isfolder(path)
        mkdir(path)
    end
    cd (path)
    
    name = ['genes_GO_term_analysis_state_', gene_state_no];
    writecell (genes, name)
    cd (currentPath)
end
end 

