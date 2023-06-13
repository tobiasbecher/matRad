function [cstIt,cstIt2,PriorityList,PriorityList2] = matRad_cstObjectiveUpdate(cstIt,cstIt2,PriorityList,PriorityList2,goal,fopt)
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
idx = PriorityList{1,1};
objective = cstIt{idx,6}{end};

if goal < fopt && fopt > 1e-4 %not met
    slack = 1.03;
    goal = fopt*slack;
    %add to cst2 as constraint
else
    PriorityList2(end+1,1:2) = PriorityList(1,:); %add to priority list in iteration 2 slack?
    PriorityList2{end,3} =  numel(cstIt2{idx,6})+1; %add index j in cstIt2{i,j} that it is stored in for second iteration
end

%always add objective as constraint for next iteration
cstIt2{idx,6}{end+1} = objective.turnIntoLexicographicConstraint(goal);

%add objective always as constraint for next iteration
constr = objective.turnIntoLexicographicConstraint(goal);
cstIt{idx,6}{end} = constr;

%Add to cst and update priority list
PriorityList(1,:) = [];

end        

