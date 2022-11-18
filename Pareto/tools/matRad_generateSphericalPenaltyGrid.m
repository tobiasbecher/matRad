function [penVal,sample] = matRad_generateSphericalPenaltyGrid(nPoints,sizes)
% matRad function that generates points for Pareto analysis, 
% distributed on a sphere based on random sampling
% call
%   [penVal,sample] = matRad_generateSphericalPenaltyGrid(nPoints,nDims)
%
% input
%   nPoints:    Number of points for penalty grid
%   sizes:      Array storing the number of objectives for each VOI
%
% output
%   sample:     Coordinates of points on plane
%   penVal:     Cell array that stores the penalties for all iterations
%
% References
%  DOI: 10.1214/009117904000000874
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


randNums= abs(normrnd(0,1,[nPoints,sum(sizes)]));
r = sqrt(sum(randNums.^2,2));
sample = randNums./r;
sample = matRad_AdjustedTravellingSalesman(sample);
sample = sample*1000;
%generate cell array for easier looping
penVal = cell(size(sizes));
k = 1;
for i=1:numel(sizes)
    VOIData = sample(:,k:k-1+sizes(i));
    %VOIData
    penVal{i} = VOIData;
    k = k+sizes(i);
end
