function cst = matRad_updatecst(cst,newPenalties,groups)
% matRad function that updates the objective penalties in the cst struct
% input
%   cst:        cst struct 
%   groups:     Cell array storing the grouping of the objectives
%   (e.g{[1,2,4],3,5}
% output
%   cst:        Modified cst struct
%   
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('groups','var')
    groups = num2cell(1:numel(newPenalties)); %check that there are as many objectives?
end
%maybe do something if groups are not specified

%check that there are as many newPenalties as groups
if numel(groups) ~= numel(newPenalties)
    matRad_cfg.dispError('Number of elements in penalty vector does not match up with groups!\n');
end


for groupidx = 1:numel(groups)

    el = groups{groupidx};
    %if the group is not represented by a cell array already, modify
    
    %loop over objectives in group
    for objNumInGroup = 1:numel(el)
        objidx = el(objNumInGroup); %How many objective is this

        %update corresponding value
        curObjCount = 0;
        for i = 1:size(cst,1) % loop over cst
            for j = 1:numel(cst{i,6})
                %check whether dose objective or constraint
                obj = cst{i,6}{j};
                if isstruct(cst{i,6}{j})
                    obj =  matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
                end 
                if contains(class(obj),'DoseObjectives')
                    curObjCount = curObjCount + 1; %update current objective number
                    if curObjCount == objidx
                       cst{i,6}{j}.penalty = newPenalties(groupidx)*100; %!penalties modified by 100
                    end
                end
            end
        end
    end
        
end
%