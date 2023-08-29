classdef matRad_PriorityList1 < matRad_PriorityClass
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

        function obj = matRad_PriorityList1() 
            obj = obj@matRad_PriorityClass();
        end

        function [cst,Priority]= modifyCst(obj,cst)
            [Priority,nextObjectives] = obj.nextPriority();
            for i = 1:numel(nextObjectives)
                cst{nextObjectives{i}.cstIdx,6}{end+1}  = nextObjectives{i}.objective;
                %set number where objective is stored in cst
                nextObjectives{i}.setVOIIdx(numel(cst{nextObjectives{i}.cstIdx,6}));
            end
        end

        function [cst1,cst2,PriorityList2] = updateStep(obj,cst1,cst2,PriorityList2,fopt)
            %Change objectives to constraints

            [Priority,LexiObjectives] = obj.nextPriority();
            %remove from current PriorityList
            obj.numOfObj = obj.numOfObj + numel(LexiObjectives);

            %change objectives in cst to constraints and add them to PriorityList2
            for i = 1:numel(LexiObjectives)
                LexiObjective = LexiObjectives{i};
                LexiObjective.achievedValue = fopt(i);
                bound  = LexiObjective.goalValue;   
                if bound < fopt(i) && fopt(i) > 1e-4 %not met
                    bound = fopt(i)*obj.slackVariable;
                    %add to cst2 as constraint
                else
                    %objectives are used in PriorityList2
                    PriorityList2.GoalList{end+1} = LexiObjective;
                    PriorityList2.Priorities = [PriorityList2.Priorities,Priority];
                end
                %update both csts 
                cst1{LexiObjective.cstIdx,6}{LexiObjective.getVOIIdx()} = LexiObjective.objective.turnIntoLexicographicConstraint(bound);
                cst2{LexiObjective.cstIdx,6}{LexiObjective.getVOIIdx()} = LexiObjective.objective.turnIntoLexicographicConstraint(bound);
            end
        end

        function [objectives,skip] = fastObjectivecalc(obj,dij,cst,optiProb,wInit)
            %fast check if objectives are already met for next iteration
            %and do not have to be recalculated
            
            [Priority,LexiObjectives] = obj.nextPriority();

            %preassigning
            bounds = zeros(numel(LexiObjectives),1);
            objidxs = [];
            objectives = {};
            for i = 1:numel(LexiObjectives)
                bounds(i) = LexiObjectives{i}.goalValue;  
                objidxs = [objidxs;LexiObjectives{i}.cstIdx,LexiObjectives{i}.getVOIIdx()];
                objectives{end+1} = LexiObjectives{i}.objective;
            end

            %modify optimizationProblem
            optiProb.objidx = objidxs;
            optiProb.objectives = objectives;
            objectives = matRad_objectiveFunctions(optiProb,wInit,dij,cst);

            skip = all(objectives < bounds);
        end

        function cst = generateOverlapcst(obj,cst)
            cst = obj.generateBasecst(cst); %add constraints
            for i = 1:numel(obj.GoalList)
                cst{obj.GoalList{i}.cstIdx,6}{end+1} = obj.GoalList{i}.objective;
            end
        end
    end
end