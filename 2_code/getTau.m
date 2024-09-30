function tau = getTau(data, rescale_data, interval_x, interval_y)

% Calculate tau for CP model fitting. 
%
% data             an n x 2 array. First column is the x-variable (time), 
%                  second column is the y-variable (gene expression). Can 
%                  be bootstrapped data.
% rescale_data     whether to rescale the data. Default false
% interval_x       an array of two elements. E.g. [0 1]
% interval_y       an array of two elements. E.g. [0 1]
%
% Chen Chen. Last update: 2024-08-04
% Rosemary Yu. Last update: 2024-08-08

if nargin < 2
    rescale_data = false;
end

if rescale_data && ~exist('interval_x', 'var') && ~exist('interval_y', 'var')
    error('rescaling interval needed for x or y variables (or both)')
end


BS2 = data;
if rescale_data && exist('interval_x', 'var')
    BS2 (:,1) = rescale (BS2 (:,1), interval_x(1), interval_x(2));
end
if rescale_data && exist('interval_y', 'var')
    BS2 (:,2) = rescale (BS2 (:,2), interval_y(1), interval_y(2));
end
%BS = BS2;

dist = [];
distp = [];
for i = 2:length (BS2)
    dist (i-1,:) = sqrt ((BS2 (i,1)-BS2(i-1,1))^2+(BS2(i,2)-BS2(i-1,2))^2);
end
dist2 = dist;
for i = 2:length (dist)
    dist2 (i,1) = dist (i,1) + dist2 (i-1,1);
end
s = sum (dist);
distp =  dist2 / s;
tau = cat (1, 1e-5, distp (1:length (distp)-1,:), 1-1e-5);


end