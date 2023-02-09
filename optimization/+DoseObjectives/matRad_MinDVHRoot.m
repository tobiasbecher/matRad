classdef matRad_MinDVHRoot < DoseObjectives.matRad_DoseObjective
% matRad_MinDVHRoot Implements a penalized MinDVH objective
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
    
    properties (Constant)
        name = 'Min DVH Root';
        parameterNames = {'d', 'V^{min}'};
        parameterTypes = {'dose','numeric'};
    end
    
    properties
        parameters = {60,95};
        penalty = 1;
        objValue = 0;
    end
    
    methods
        function obj = matRad_MinDVHRoot(penalty,dRef,vMinPercent)
            
            % if we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DoseObjectives.matRad_DoseObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin >= 3 && isscalar(vMinPercent)
                    obj.parameters{2} = vMinPercent;
                end
                
                if nargin >= 2 && isscalar(dRef)
                    obj.parameters{1} = dRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
            
        end        
        %% Calculates the Objective Function value
        function fDose = computeDoseObjectiveFunction(obj,dose)                       
            % get reference Volume
            refVol = obj.parameters{2}/100;
            
            % calc deviation
            deviation = dose - obj.parameters{1};

            % calc d_ref2: V(d_ref2) = refVol
            d_ref2 = matRad_calcInversDVH(refVol,dose);

            
            deviation(dose > obj.parameters{1} | dose < d_ref2) = 0;
   
            % claculate objective function
            fDose = sqrt(1/numel(dose) * (deviation'*deviation));
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad   = computeDoseObjectiveGradient(obj,dose)
            % get reference Volume
            refVol = obj.parameters{2}/100;
            
            % calc deviation
            deviation = dose - obj.parameters{1};
            
            % calc d_ref2: V(d_ref2) = refVol
            d_ref2 = matRad_calcInversDVH(refVol,dose);
            
            deviation(dose > obj.parameters{1} | dose < d_ref2) = 0;

            % calculate delta
            fDoseGrad = 1/sqrt(numel(dose)*(deviation'*deviation)) * deviation;
        end
    end
    
end