function [PriorityList,cstConstr] = matRad_PriorityList(cst)
% Given a cst file return a priority list of the objectives as well as a
% list of constraints with corresponding indices
% 
% call
%
% input
%   cst:            modified matRad cst struct. 
%
% output
%   PriorityList:   cell array storing the objectives and their respective
%   indices in the cst for all objectives in order of priority
%   cstConstr:      Modified cst file that only contains constraints
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
%preallocation
Priorities = [];
PriorityList = {};
ConstraintList = {};
PriorityCounter = 1;
ConstraintCounter = 1;
for i = 1:size(cst,1) % loop over cst
    
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
    
        for j = 1:numel(cst{i,6})
            %check whether dose objective or constraint
            obj = cst{i,6}{j};
            if isstruct(cst{i,6}{j})
                obj =  matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
            end
            if contains(class(obj),'DoseObjectives')
                %add objectives to priority list
                %does not handle same priorities right now
                PriorityList{PriorityCounter,1} = i;
                PriorityList{PriorityCounter,2} = obj;

                Priorities = [Priorities,obj.penalty];
                PriorityCounter = PriorityCounter + 1;

            else %constraint
                %add constraints to respective list
                ConstraintList{ConstraintCounter,1} = i;
                ConstraintList{ConstraintCounter,2} = obj;
                ConstraintCounter = ConstraintCounter +1;
            end

        end
    end
end
%sort PriorityList based on priorities
[Priorities,I] = sort(Priorities);
PriorityList = PriorityList(I,:);


%modify cst to only store initial hard constraints
cst(:,6) = cell(1);
for i = 1:size(ConstraintList,1)
    cst{ConstraintList{i,1},6}{end+1} = ConstraintList{i,2};
end
cstConstr = cst;
end