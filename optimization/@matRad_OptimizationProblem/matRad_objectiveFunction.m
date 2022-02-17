function f = matRad_objectiveFunction(optiProb,w,dij,cst)
% matRad IPOPT objective function wrapper
% 
% call
%   f = matRad_objectiveFuncWrapper(optiProb,w,dij,cst)
%
% input
%   optiProb: matRad optimization problem
%   w:        beamlet/ pencil beam weight vector
%   dij:      matRad dose influence struct
%   cst:      matRad cst struct
%
% output
%   f: objective function value
%
% References
%   -
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

% get current dose / effect / RBExDose vector
%d = optiProb.matRad_backProjection(w,dij);
optiProb.BP = optiProb.BP.compute(dij,w);
d = optiProb.BP.GetResult();

if ~isempty(optiProb.BP_LET)
    optiProb.BP_LET = optiProb.BP_LET.compute(dij,w);
    LET = optiProb.BP_LET.GetResult();
end

if ~isempty(optiProb.BP_DADRfixed)
    optiProb.BP_DADRfixed = optiProb.BP_DADRfixed.compute(dij,w);
    DADR = optiProb.BP_DADRfixed.GetResult();
end

if ~isempty(optiProb.BP_XBDDADR)
    optiProb.BP_XBDDADR = optiProb.BP_XBDDADR.compute(dij,w);
    XBDDADR = optiProb.BP_XBDDADR.GetResult();
end

if ~isempty(optiProb.BP_XBDLET)
    optiProb.BP_XBDLET = optiProb.BP_XBDLET.compute(dij,w);
    XBDLET = optiProb.BP_XBDLET.GetResult();
end


% Initialize f
f = 0;

% compute objectiveective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            
            objective = cst{i,6}{j};
            
            % only perform gradient computations for objectiveectives
            %if isempty(strfind(objective.type,'constraint'))
            if isa(objective,'DoseObjectives.matRad_DoseObjective')

                % if we have effect optimization, temporarily replace doses with effect
                if (~isequal(objective.name, 'Mean Dose') && ~isequal(objective.name, 'EUD')) &&...
                    (isa(optiProb.BP,'matRad_EffectProjection') && ~isa(optiProb.BP,'matRad_VariableRBEProjection')) 
                    
                    doses = objective.getDoseParameters();
                
                    effect = cst{i,5}.alphaX*doses + cst{i,5}.betaX*doses.^2;
                    
                    objective = objective.setDoseParameters(effect);
                end
                
                % if conventional opt: just sum objectiveectives of nominal dose
                %if strcmp(cst{i,6}{j}.robustness,'none')

                    d_i = d{1}(cst{i,4}{1});

                    f = f + objective.computeDoseObjectiveFunction(d_i);
                    
                %end
            
            end
            
            % only perform gradient computations for objectiveectives
            %if isempty(strfind(objective.type,'constraint'))
            if isa(objective,'LETObjectives.matRad_LETObjective')
              
                % if conventional opt: just sum objectiveectives of nominal LET
                %if strcmp(cst{i,6}{j}.robustness,'none')

                    LET_i = LET{1}(cst{i,4}{1});

                    f = f + objective.computeLETObjectiveFunction(LET_i);
                    
                %end
            
            end

            if isa(objective,'DADRObjectives.matRad_DADRObjective')
              
                % if conventional opt: just sum objectiveectives of nominal LET
                %if strcmp(cst{i,6}{j}.robustness,'none')

                    DADR_i = DADR{1}(cst{i,4}{1});

                    f = f + objective.computeDADRObjectiveFunction(DADR_i);
                    
                %end
            
            end

            if isa(objective,'XBDDADRObjectives.matRad_XBDDADRObjective')
              
                % if conventional opt: just sum objectiveectives of nominal LET
                %if strcmp(cst{i,6}{j}.robustness,'none')

                    XBDDADR_i = XBDDADR{1}(cst{i,4}{1});

                    f = f + objective.computeXBDDADRObjectiveFunction(XBDDADR_i);
                    
                %end
            
            end

            if isa(objective,'XBDLETObjectives.matRad_XBDLETObjective')
              
                % if conventional opt: just sum objectiveectives of nominal LET
                %if strcmp(cst{i,6}{j}.robustness,'none')

                    XBDLET_i = XBDLET{1}(cst{i,4}{1});

                    f = f + objective.computeXBDLETObjectiveFunction(XBDLET_i);
                    
                %end
            
            end
       
        end
            
    end
    
end
