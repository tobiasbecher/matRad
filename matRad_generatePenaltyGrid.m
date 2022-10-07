function penalties = matRad_generatePenaltyGrid(pen)
% matRad function to create an easily loopable penalty grid
% 
%
% input
%   pen:        STRUCT WITH DIFFERENT PENALTIES FOR PARETOSURFACE
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
%create cell array for 
grid = cell(size(pen));
[grid{:}] = ndgrid(pen{:});
%reshape the grid to be looped over
penalties = zeros(size(grid{1}(:),1),size(grid,2));
for i=1:size(penalties,2)
    penalties(:,i) = grid{i}(:);
end