function [cstIt,cstIt2,PriorityList,idx,goal] = matRad_UpdateCst(cstIt,PriorityList,idx,fopt)
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
%Turn last used objective into constraint
objective = cstIt{idx,6}{end};
constr = objective.turnIntoLexicographicConstraint(fopt);
cstIt{idx,6}{end} = constr;

%Add next objective from priority list
nObj = PriorityList{1,2};
[nObj,goal] = nObj.SetAsLexicographic();

cstIt{PriorityList{1,1},6}{end+1} = nObj;
idx = PriorityList{1,1};
PriorityList(1,:) = [];
end