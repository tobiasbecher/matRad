function [cstIt,PriorityList,goal] = matRad_cstNextObjective(cstIt,PriorityList)
% Given a cst file containing constraints add the next objective from the
% priority list and update the PriorityList
% 
% call
%
% input
%   cst:            modified matRad cst struct. 
%   PriorityList:   PriorityList cell array.
%
% output
%   priorities: Array storing the priorities of the objectives
%   idxs:       Array storing the corresponding cst indices
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%modify objective for lexicographic optimization
objective = PriorityList{1,2};
[objective,goal] = objective.SetAsLexicographic();

%Add to cst and update priority list
cstIt{PriorityList{1,1},6}{end+1} = objective;
end

