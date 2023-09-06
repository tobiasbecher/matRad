function resultGUIVis = matRad_convert2pPlans(resultGUIs,pln)
% Converts the resultGUI for a mixed Mod plan to a resultGUI that can be
% visualized by the GUI
%
% call
%   resultGUI = matrad_convert2pPlans(resultGUI,pln)
%
% input
%   resultGUIs:  Cell array storing all "plans" of a given step
%   pln:        struct storing the mixed Mod plan
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

for i = 1:numel(resultGUIs)
    resultGUI = resultGUIs{i};
    istr = num2str(i);
    resultGUIVis.(['physicalDose_',istr]) = zeros(size(resultGUI{1}.physicalDose)); 

    for j = 1:numel(pln.originalPlans)
        phDose = resultGUI{j}.physicalDose*pln.originalPlans(j).numOfFractions;
        resultGUIVis.(['physicalDose_',istr]) = resultGUIVis.(['physicalDose_',istr]) + phDose;
        resultGUIVis.(['physicalDose_',istr,'_', pln.originalPlans(j).radiationMode]) = phDose;
    end
end