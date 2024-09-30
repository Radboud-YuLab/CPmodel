function [nCP,CPloc] = writeOptCP(out_opt, savedir)
% saves the output in two .csv files for user friendliness. 
%
% nCP           m x 5 cell array with the following:
%               col 1: gene name
%               col 2: number of CPs
%               col 3: RSS
%               col 4: nCP-penalized RSS
%               col 5: R squared
% CPloc         n x 5 cell array with the following:
%               col 1: gene name
%               col 2: CP number
%               col 3: time coordinate of that CP
%               col 4: expression coordinate of the CP
%               col 5: state (initial, intermediate or final)
% out_opt       a j x 6 cell array collecting the second output of the
%               getOptCP function for m number of genes. Should contain
%               the following:
%               col 1: gene name
%               col 2: number of CPs
%               col 3: array of CPs. Each row is one CP in (x,y)
%               coordinates.
%               col 4: RSS
%               col 5: nCP-penalized RSS
%               col 6: R squared
% savedir       directory to save the files
%
% Chen Chen. Last updated: 2024-09-06
nCP = {};
CPloc = {};
for i = 1:length (out_opt);
    nCP {i,1} = out_opt {i,1};
    nCP {i,2} = out_opt {i,2};
    nCP {i,3} = out_opt {i,4};
    nCP {i,4} = out_opt {i,5};
    nCP {i,5} = out_opt {i,6};
end
tnCP=cell2table(nCP,'VariableNames',{'Gene name','number of CPs','RSS','nCP-penalized RSS','R squared'});

for i = 1:length (out_opt);
    subCPloc = {};
    j = [];
    if isempty (out_opt {i,2})
        subCPloc {1,1} = out_opt {i,1};
        for o = 2:5
        subCPloc {1,o} = 'N/A';
        end
    else
        ncp = out_opt {i,2};
        for j = 1:ncp
            subCPloc {j,1} = out_opt {i,1};
            subCPloc {j,2} = j;
            subCPloc {j,3} = out_opt{i, 3}(j,1);
            subCPloc {j,4} = out_opt{i, 3}(j,2);
            subCPloc {j,5} = 'intermediate';
        end
        subCPloc {1,5} = 'initial';
        subCPloc {height (subCPloc), 5} = 'final';
    end
    CPloc = cat (1, CPloc, subCPloc);
end
tCPloc=cell2table(CPloc,'VariableNames',{'Gene name','CP number','CP_t','CP_exp','state'});

if exist('savedir', 'var')
    %make directory
    currentPath = pwd;
    if ~isfolder(savedir)
        mkdir(savedir)
    end
    cd (savedir)

    writetable (tnCP, 'nCP.csv');
    writetable (tCPloc, 'CPloc.csv');
    cd (currentPath)
else
    display ('No save directory was given')
    cd (currentPath)
end

end

