function [penVal,penGrid] = matRad_generateComplexPenaltyGrid(VOI, penalties)
% matRad function to create an easily loopable penalty grid
% 
%
% input
%   VOI 
%
% output
%   penalty:    A grid to loop over; gridToLoop(i,:) to access the current
%               penalty values
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%VOI = penaltyMap.keys; %!!! Sorts keys alphabetically
%penalties = penaltyMap.values;

%create array storing number of objectives/penalties for each VOI
sizes = zeros(1,numel(VOI));
for i=1:numel(VOI)    
    indPen = penalties{i};
    sizes(i) = size(indPen,2);
end
sizes
% convert grid with VOIs with possibly more than one objective function
singPen = cell(1,sum(sizes));
k = 1;
for i = 1:numel(penalties)
    penalty = penalties{i};
    for j = 1:numel(penalty)
        singPen{k} = cell2mat(penalty{j});
        k = k + 1;
    end
end
%use simple grid function to generate Penalty Grid
penGrid = matRad_generatePenaltyGrid(singPen);

%renormalizing penalty values so all rows add up to 1000
penGrid = bsxfun(@rdivide, penGrid, (sum(penGrid,2))/1000);

%remove rows that only contain nan values (happens when all penalties 
%are 0 (might be possible to simply remove the first row after unique)
penGrid = penGrid(all(~isnan(penGrid),2),:);

%eliminating same rows
penGrid = unique(penGrid,'rows','stable');

%reshape Penalty Grid to respect multiobjectiv VOIs
penVal = cell(size(VOI));
k = 1;
for i=1:numel(sizes)
    VOIData = penGrid(:,k:k-1+sizes(i));
    penVal{i} = VOIData;
    k = k+sizes(i);
end