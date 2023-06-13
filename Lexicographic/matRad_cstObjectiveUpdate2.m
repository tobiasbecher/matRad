function [cstIt2,PriorityList2] = matRad_cstObjectiveUpdate2(cstIt2,PriorityList2,fopt)
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
idx = PriorityList2{1,1};
idx2 = PriorityList2{1,3};
objective = cstIt2{idx,6}{idx2};
%{
if goal < fopt && fopt > 1e-4 %not met
    slack = 1.03;
    goal = fopt*slack;
end
%}

%add objective always as constraint for next iteration
slack = 1.03;
constr = objective.turnIntoLexicographicConstraint(fopt*slack);
cstIt2{idx,6}{idx2} = constr;

%Add to cst and update priority list
PriorityList2(1,:) = [];

end        

