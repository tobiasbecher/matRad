function orderedPenPoints = matRad_orderPoints(penPoints)
% matRad function to reorder penaltyGrid to be more favorable for
% optimization
%
% input
%   penPoints:          matrix containing the penalty Points to be reordered
%
% output
%   orderedpenPoints:   reordered matrix
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
orderedPenPoints = zeros(size(penPoints));
refPoint = zeros(1,size(penPoints,2));
refPoint(1,1) = 1;
for i = 1:size(penPoints,1)
    dist = sum((penPoints-refPoint).^2,2);
    [M,I] =  min(dist);
    orderedPenPoints(i,:) = penPoints(I,:);
    refPoint = penPoints(I,:);
    penPoints(I,:) = [];
end