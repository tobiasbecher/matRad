function [resultGUI,resultGUIs,resultGUIs2,cstIt] =  matRad_2pecOptimization(dij,cst,pln,wInit)
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
    slack = 1.03;
    %% First Phase: Find aspiration levels



    [PriorityList,cstIt] = matRad_PriorityList(cst); %extract the priorities and indices from the objective
    cstIt2 = cstIt; %copy for second iteration
    PriorityList2 = {};


    %add first objective to cst and remove it from PriorityList
    [cstIt,PriorityList,goal] = matRad_cstNextObjective(cstIt,PriorityList);
    goal = goal/pln.numOfFractions; % works for common objectives, not sure about dvh

    %iterative optimization
    
    for i = 1:size(PriorityList,1)-1
        if exist('wInit','var')
            resultGUI = matRad_fluenceOptimization(dij,cstIt,pln,wInit);
        else
            resultGUI = matRad_fluenceOptimization(dij,cstIt,pln);
        end

        %change objective in cstIt to constraint and add it to 2nd
        %constraint or priority list depending on outcome
        [cstIt,cstIt2,PriorityList,PriorityList2] = matRad_cstObjectiveUpdate(cstIt,cstIt2,PriorityList,PriorityList2,goal,resultGUI.fopt);

        %add next objective from priority list
        [cstIt,PriorityList,goal] = matRad_cstNextObjective(cstIt,PriorityList);

        wInit = resultGUI.w;
        resultGUIs{end+1} = resultGUI;
        goal = goal/pln.numOfFractions; %should not be universal (e.g for sq objs)
    end
    
        %last run of optimization
    if exist('wInit','var')
        resultGUI = matRad_fluenceOptimization(dij,cstIt,pln,wInit);
    else
        resultGUI = matRad_fluenceOptimization(dij,cstIt,pln);
    end
    %final updates
    wInit = resultGUI.w;
    [cstIt,cstIt2,PriorityList,PriorityList2] = matRad_cstObjectiveUpdate(cstIt,cstIt2,PriorityList,PriorityList2,goal,resultGUI.fopt);
    
    resultGUIs{end+1} = resultGUI;

    %% Second phase: 
    resultGUIs2 = {};
    [cstIt2,PriorityList2,goal] = matRad_cstNextObjective2(cstIt2,PriorityList2);

    for i = 1:size(PriorityList2,1)-1
        if exist('wInit','var')
            resultGUI = matRad_fluenceOptimization(dij,cstIt2,pln,wInit);
        else
            resultGUI = matRad_fluenceOptimization(dij,cstIt2,pln);
        end
    
        %change objective in cstIt to constraint and add it to 2nd
        %constraint or priority list depending on outcome
        [cstIt2,PriorityList2] = matRad_cstObjectiveUpdate2(cstIt2,PriorityList2,resultGUI.fopt);
    
        %add next objective from priority list
        [cstIt2,PriorityList2,goal] = matRad_cstNextObjective2(cstIt2,PriorityList2);
    
        wInit = resultGUI.w;
        resultGUIs2{end+1} = resultGUI;
        goal = goal/pln.numOfFractions; %should not be universal (e.g for sq objs)
    end
    
    %last run, so wInit should always be fdefined

    %[cstIt2,PriorityList2,goal] = matRad_cstNextObjective(cstIt2,PriorityList2);

    resultGUI = matRad_fluenceOptimization(dij,cstIt2,pln,wInit);
    resultGUIs2{end+1} = resultGUI;
end