function [resultGUIs,resultGUIs2,cst1,cst2,PriorityList2]=  matRad_2pecOptimizationNew(PriorityList,dij,ct,cst,pln,wInit)
    % Lexicographic optimization
    % 
    % call
    %
    % input
    %   dij:        matRad dij struct
    %   cst:        modified matRad cst struct. 
    %   pln:        matRad pln struct
    %   wInit:      (optional) custom weights to initialize problems
    %
    % output
    %
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
    resultGUIs = {};
    resultGUIs2 = {};
    PriorityList2 = matRad_PriorityList2();
    PriorityList2.ConstraintList = PriorityList.ConstraintList; %same constraints
    
    %initial update for cst with all objectives and constraints for overlap
    cst = PriorityList.generateOverlapcst(cst);
    cst = matRad_setOverlapPriorities(cst);
    cst1 = PriorityList.generateBasecst(cst);
    cst2 = cst1;

    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.propOpt.defaultAccChangeTol = 1e-6;
    %add first objective(s) to cst
    cst1 =  PriorityList.modifyCst(cst1);
    while PriorityList.numOfObj <= numel(PriorityList.GoalList)
        if exist('wInit','var')
            %check if objective can be skipped
            cstR = matRad_resizeCstToGrid(cst1,dij.ctGrid.x,  dij.ctGrid.y,  dij.ctGrid.z,...
                                 dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
            [objectives,FastCalc] = PriorityList.fastObjectivecalc(dij,cstR,optiProb,wInit);

            if ~FastCalc
                [resultGUI,optimizer,optiProb] = matRad_fluenceOptimizationJO(dij,cst1,pln,wInit); %should be choseable
                objectives = resultGUI{end}.objectives;
            else
                resultGUI{end}.objectives = objectives;
            end
            
        else
            [resultGUI,optimizer,optiProb]  = matRad_fluenceOptimizationJO(dij,cst1,pln);    
            objectives = resultGUI{end}.objectives;
        end
        wInit = resultGUI{1}.w;
        for i= 2:numel(resultGUI)-1
            wInit = [wInit;resultGUI{i}.w];
        end
        
        resultGUIs{end+1} = resultGUI;
        % exchange previously optimized objectives to constraints and remove from priorityList
        [cst1,cst2,PriorityList2] = PriorityList.updateStep(cst1,cst2,PriorityList2,objectives); 
        
        %add next objective

        if PriorityList.numOfObj <= numel(PriorityList.GoalList)
            cst1 = PriorityList.modifyCst(cst1);
        end
    end
    %step2
    clear wInit;
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.propOpt.defaultAccChangeTol = 1e-6;
    cst2 = PriorityList2.modifyCst(cst2); %should set objective in appropriate spot -> use VOIIdx
    while PriorityList2.numOfObj <= numel(PriorityList2.GoalList)
        if exist('wInit','var')
            resultGUI = matRad_fluenceOptimizationJO(dij,cst2,pln,wInit); %should be choseable
        else
            resultGUI = matRad_fluenceOptimizationJO(dij,cst2,pln);
        end

        wInit = resultGUI{1}.w;
        for i= 2:numel(resultGUI)-1
            wInit = [wInit;resultGUI{i}.w];
        end
    

        resultGUIs2{end+1} = resultGUI;
        [cst2] = PriorityList2.updateStep(cst2,resultGUI{end}.objectives); 

        if PriorityList2.numOfObj <= numel(PriorityList2.GoalList)
            cst2 = PriorityList2.modifyCst(cst2);
        end

    end
end