classdef(Abstract) matRad_PriorityClass < handle
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
    properties
        Priorities;
        ConstraintList;
        GoalList;
        slackVariable = 1.03;
        numOfObj;
    end

        
    methods
        function obj = matRad_PriorityClass()
            %add optional slack Variable adjustment
            obj.Priorities = [];
            obj.GoalList = {};
            obj.ConstraintList = {};
            obj.numOfObj = 1;

        end
        
        function addObjective(obj,Priority,objective,goal,cstIdx)
            obj.Priorities = [obj.Priorities,Priority];
            obj.GoalList{end+1} = matRad_PriorityListObjective(objective,goal,cstIdx);
            %sort by priority
            [obj.Priorities,I] = sort(obj.Priorities);
            obj.GoalList = obj.GoalList(I);
        end

        function addConstraint(obj,constraint,cstIdx)
            obj.ConstraintList{end+1} = matRad_PriorityListConstraint(constraint,cstIdx);
        end


        function cst = generateBasecst(obj,cst)
            %generate a bare bone cst struct for prioritized optimization containing only the 
            % hard constraints
            cst(:,6) = cell(1);
            for i = 1:numel(obj.ConstraintList)
                cst{obj.ConstraintList{i}.cstIdx,6}{end+1} = obj.ConstraintList{i}.constraint;
            end
        end

        function [Priority,LexiObjectives] = nextPriority(obj) %get next element in 1st prioList
            Priority = obj.Priorities(obj.numOfObj); %get next highest priority and remove all objectives from list
            LexiObjectives = {};
            i = obj.numOfObj;
            while i <= numel(obj.Priorities) && obj.Priorities(i) == Priority 
                LexiObjectives{end+1} = obj.GoalList{i};
                i = i + 1;
            end
        end

        function printPriorityList(obj,cst)
            %create data
            goals = [];
            VOINames= {};
            Objectives = {};
            AchievedValues = [];
            AchievedValues2 = [];
            for i = 1:numel(obj.GoalList)
                idx = obj.GoalList{i}.cstIdx;
                goals = [goals;obj.GoalList{i}.goalValue];
                VOINames = [VOINames;cst{idx,2}];
                Objectives = [Objectives;obj.GoalList{i}.objective.name];
                AchievedValues = [AchievedValues;obj.GoalList{i}.achievedValue];
                AchievedValues2 = [AchievedValues2;obj.GoalList{i}.achievedValue2];
            end
            Priority = obj.Priorities';
            T = table(Priority,VOINames,Objectives,goals,AchievedValues,AchievedValues2);
            fig = uifigure();
            uitable(fig,"Data",T);
        end
            

    end

    methods (Abstract)
        
        modifyCst
        updateStep  
    end

end