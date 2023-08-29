classdef matRad_PriorityList2 < matRad_PriorityClass
% matRad_MeanDose Implements a penalized MeanDose objective
%   See matRad_DoseObjective for interface description
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods 

        function obj = matRad_PriorityList2() 
            obj = obj@matRad_PriorityClass();
        end

        function [cst,Priority]= modifyCst(obj,cst)
            [Priority,nextObjectives] = obj.nextPriority();
            for i = 1:numel(nextObjectives)
                %use previously allocated position
                cst{nextObjectives{i}.cstIdx,6}{nextObjectives{i}.getVOIIdx} = nextObjectives{i}.objective;
            end

        end

        function [cst2] = updateStep(obj,cst2,fopt)
            %change the objective to constraint
            [Priority,LexiObjectives] = obj.nextPriority();
            %remove objective(s) from current PriorityList

            obj.numOfObj = obj.numOfObj + numel(LexiObjectives);

            for i = 1:numel(LexiObjectives)
                LexiObjective.achievedValue = fopt(i);
                LexiObjective = LexiObjectives{i};
                bound = fopt(i)*obj.slackVariable;
                cst2{LexiObjective.cstIdx,6}{LexiObjective.getVOIIdx()} = LexiObjective.objective.turnIntoLexicographicConstraint(bound);
            end
        end
    end
end