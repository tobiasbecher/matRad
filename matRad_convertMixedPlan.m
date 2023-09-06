function resultGUIVis = matRad_convertMixedPlan(resultGUI,pln,bound)
% Converts the resultGUI for a mixed Mod plan to a resultGUI that can be
% visualized by the GUI
%
% call
%   resultGUI = matrad_convertMixedPlan(resultGUI,pln)
%   resultGUI = matrad_convertMixedPlan(resultGUI,pln,bound)
%
% input
%   resultGUI:  struct storing the results of the mixed Mod optimization
%   pln:        struct storing the mixed Mod plan
%   bound:      (optional) set all values above bound to bound value
%
% output
%   resultGUI: matRad result struct
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

    if nargin < 3
        bound = 1e9;
    end
        
    %resultGUIVis = struct();
    %% Store initial fields

    for i = 1:numel(pln.originalPlans)
        fn = fieldnames(resultGUI{i});
        for j = 1:numel(fn)
            resultGUIVis.([fn{j},'_',pln.originalPlans(i).radiationMode]) = resultGUI{i}.([fn{j}]);
        end
    end
    
    %% Store combined fields
    
    %physical Dose
    resultGUIVis.physicalDose = zeros(size(resultGUI{1}.physicalDose)); 
    if isfield(resultGUI{1},'RBExD')
        resultGUIVis.RBExD = zeros(size(resultGUI{1}.physicalDose));
    end



    for i = 1:numel(pln.originalPlans)
        phDose = resultGUI{i}.physicalDose*pln.originalPlans(i).numOfFractions;
        phDose(phDose> bound) = bound;
        resultGUIVis.physicalDose = resultGUIVis.physicalDose + phDose;

        if isfield(resultGUI{1},'RBExD')
            RBE = resultGUI{i}.RBExD*pln.originalPlans(i).numOfFractions;
            RBE(RBE > bound) = bound;
            resultGUIVis.RBExD = resultGUIVis.RBExD + RBE;
        end
    end


    


    %calculate relative values
    if numel(resultGUI) > 2
        for i = 1:numel(pln.originalPlans)
            phDose = resultGUI{i}.physicalDose*pln.originalPlans(i).numOfFractions;
            phDose(phDose> bound) = bound;
            resultGUIVis.(['relativeDose', pln.originalPlans(i).radiationMode]) = phDose./resultGUIVis.physicalDose;
        end
    
        resultGUIVis.FracDose =  resultGUI{2}.physicalDose*pln.originalPlans(2).numOfFractions;
        resultGUIVis.FracDose =  resultGUIVis.FracDose./ resultGUI{1}.physicalDose*pln.originalPlans(1).numOfFractions;
        resultGUIVis.FracDose(resultGUIVis.FracDose >1e9) = 1e9;
    end
end