function addWeights = matRad_convexWeightExpansion(fVals,weights)
% matRad inverse pareto planning wrapper function
%
% input
%   fVals:      Array storing the individual objective function values
%   weights:    Array storing the beamlet weights associated to fVals
%   dij:        matRad dij construct
% output
%   addWeights: Array storing additional weights
%
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
k = convhulln(fVals);
addWeights = zeros(size(k),size(weights,2)); 
%loop over all surfaces
for i = 1:size(k,1)
    psweights = weights(k(i,:),:);
    addWeights(i,:) = mean(psweights,1); 
end
