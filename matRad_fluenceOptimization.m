function [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln,wInit)
% matRad inverse planning wrapper function
% 
% call
%   [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln)
%   [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln,wInit)
%
% input
%   dij:        matRad dij struct
%   cst:        matRad cst struct
%   pln:        matRad pln struct
%   wInit:      (optional) custom weights to initialize problems
%
% output
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%   optimizer:  Used Optimizer Object
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

matRad_cfg = MatRad_Config.instance();

% issue warning if biological optimization impossible
if sum(strcmp(pln.propOpt.bioOptimization,{'LEMIV_effect','LEMIV_RBExD'}))>0 && (~isfield(dij,'mAlphaDose') || ~isfield(dij,'mSqrtBetaDose')) && strcmp(pln.radiationMode,'carbon')
    warndlg('Alpha and beta matrices for effect based and RBE optimization not available - physical optimization is carried out instead.');
    pln.propOpt.bioOptimization = 'none';
end

% consider VOI priorities
cst  = matRad_setOverlapPriorities(cst);

haveDoseObjectives = false;
haveLETobjectives = false;
haveDADRobjectives = false;

% check & adjust objectives and constraints internally for fractionation 
for i = 1:size(cst,1)
    %Compatibility Layer for old objective format
    if isstruct(cst{i,6})
        cst{i,6} = arrayfun(@matRad_DoseOptimizationFunction.convertOldOptimizationStruct,cst{i,6},'UniformOutput',false);
    end
    for j = 1:numel(cst{i,6})
        
        obj = cst{i,6}{j};              
  
        
        %In case it is a default saved struct, convert to object
        %Also intrinsically checks that we have a valid optimization
        %objective or constraint function in the end
        if isstruct(obj)
            if strncmp(obj.className,'DoseObjective',13)
                try
                    obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
                    haveDoseObjectives = true;
                catch
                    matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!',i,j);
                end           
            
            
             elseif strncmp(obj.className,'LETObjective',12)
                try
                    obj = matRad_LETOptimizationFunction.createInstanceFromStruct(obj);
                    haveLETobjectives = true;
                catch
                    matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!',i,j);
                end
            elseif strncmp(obj.className,'DADRObjective',13)
                try
                    obj = matRad_DADROptimizationFunction.createInstanceFromStruct(obj);
                    haveDADRobjectives = true;
                catch
                    matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!',i,j);
                end
            end
        end
        
        if isa(obj,'matRad_DoseOptimizationFunction')
            obj = obj.setDoseParameters(obj.getDoseParameters()/pln.numOfFractions);
        elseif isa(obj,'matRad_LETOptimizationFunction') 
            obj = obj.setLETdParameters(obj.getLETdParameters()/pln.numOfFractions); %If it is an LET*d parameter, we need to scale it with fractions
        elseif isa(obj,'matRad_DADROptimizationFunction')
            obj = obj.setDoseParameters(obj.getDoseParameters()/pln.numOfFractions);            
        else
            matRad_cfg.dispError('Unknown objective of type %s!',class(obj));
        end
        cst{i,6}{j} = obj;        
    end
end

% Quantity availability check
if haveLETobjectives && ~isfield(dij,'mLETDose')
    matRad_cfg.dispError('LET objectives set, but no LET available in dij!');
end

if haveDADRobjectives && ~isfield(dij,'fixedCurrent')
    dij.fixedCurrent = 300; %nA
    matRad_cfg.dispWarning('DADR objectives set, but no current given in dij.fixedCurrent. Using current of %f nA!',dij.fixedCurrent);
end

% resizing cst to dose cube resolution 
cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
                                 dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);                             
                                                                                     
% find target indices and described dose(s) for weight vector
% initialization
V          = [];
doseTarget = [];
ixTarget   = [];

for i = 1:size(cst,1)
    for j = 1:size(cst{i,6},2)
        if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6}{j}) && isa(cst{i,6}{j},'matRad_DoseOptimizationFunction')
            V = [V;cst{i,4}{1}];

            %Iterate through objectives/constraints
            fDoses = [];
            for fObjCell = cst{i,6}{j}
                dParams = fObjCell.getDoseParameters();
                %Don't care for Inf constraints
                dParams = dParams(isfinite(dParams));
                %Add do dose list
                fDoses = [fDoses dParams];
            end  

            for fObjCell = cst{i,6}{j}
                dParams = fObjCell.getDoseParameters();
                %Don't care for Inf constraints
                dParams = dParams(isfinite(dParams));
                %Add do dose list
                fDoses = [fDoses dParams];
            end 

            doseTarget = [doseTarget fDoses];
            ixTarget   = [ixTarget i*ones(1,length(fDoses))];
        
        elseif isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6}{j}) && isa(cst{i,6}{j},'matRad_LETOptimizationFunction')
            V = [V;cst{i,4}{1}];
            %Iterate through objectives/constraints
            fLETs = [];
            for fObjCell = cst{i,6}{j}
                LETParams = fObjCell.getLETdParameters();
                %Don't care for Inf constraints
                LETParams = LETParams(isfinite(LETParams));
                %Add do dose list
                fLETs = [fLETs LETParams];
            end  

            for fObjCell = cst{i,6}{j}
                LETParams = fObjCell.getLETdParameters();
                %Don't care for Inf constraints
                LETParams = LETParams(isfinite(LETParams));
                %Add do dose list
                fLETs = [fLETs LETParams];
            end 

            doseTarget = [doseTarget 1];
            ixTarget   = [ixTarget i*ones(1,length(fLETs))];        
        end
    end
end
[doseTarget,i] = max(doseTarget);
ixTarget       = ixTarget(i);
wOnes          = ones(dij.totalNumOfBixels,1);

                                                         
% modified settings for photon dao
if pln.propOpt.runDAO && strcmp(pln.radiationMode,'photons')
%    options.ipopt.max_iter = 50;
%    options.ipopt.acceptable_obj_change_tol     = 7e-3; % (Acc6), Solved To Acceptable Level if (Acc1),...,(Acc6) fullfiled

end
% calculate initial beam intensities wInit
if exist('wInit','var')
    %do Nothing
elseif  strcmp(pln.propOpt.bioOptimization,'const_RBExD') && strcmp(pln.radiationMode,'protons')
    
    % check if a constant RBE is defined - if not use 1.1
    if ~isfield(dij,'RBE')
        dij.RBE = 1.1;
    end
    bixelWeight =  (doseTarget)/(dij.RBE * mean(dij.physicalDose{1}(V,:)*wOnes)); 
%    bixelWeight =  (LETTarget)/(dij.RBE * mean(dij.physicalDose{1}(V,:)*wOnes)); 
    wInit       = wOnes * bixelWeight;
        
elseif (strcmp(pln.propOpt.bioOptimization,'LEMIV_effect') || strcmp(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
                                && strcmp(pln.radiationMode,'carbon')
                            
    % retrieve photon LQM parameter
    [ax,bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,1);

    if ~isequal(dij.ax(dij.ax~=0),ax(dij.ax~=0)) || ...
       ~isequal(dij.bx(dij.bx~=0),bx(dij.bx~=0))
         matRad_cfg.dispError('Inconsistent biological parameter - please recalculate dose influence matrix!\n');
    end

    for i = 1:size(cst,1)
        
        for j = 1:size(cst{i,6},2)
            % check if prescribed doses are in a valid domain
            if any(cst{i,6}{j}.getDoseParameters() > 5) && isequal(cst{i,3},'TARGET')
                matRad_cfg.dispError('Reference dose > 10 Gy[RBE] for target. Biological optimization outside the valid domain of the base data. Reduce dose prescription or use more fractions.\n');
            end
            
        end
    end
    
    dij.ixDose  = dij.bx~=0; 
        
    if isequal(pln.propOpt.bioOptimization,'LEMIV_effect')
        
           effectTarget = cst{ixTarget,5}.alphaX * doseTarget + cst{ixTarget,5}.betaX * doseTarget^2;
           p            = (sum(dij.mAlphaDose{1}(V,:)*wOnes)) / (sum((dij.mSqrtBetaDose{1}(V,:) * wOnes).^2));
           q            = -(effectTarget * length(V)) / (sum((dij.mSqrtBetaDose{1}(V,:) * wOnes).^2));
           wInit        = -(p/2) + sqrt((p^2)/4 -q) * wOnes;

    elseif isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')
        
           %pre-calculations
           dij.gamma              = zeros(dij.doseGrid.numOfVoxels,1);   
           dij.gamma(dij.ixDose) = dij.ax(dij.ixDose)./(2*dij.bx(dij.ixDose)); 
            
           % calculate current in target
           CurrEffectTarget = (dij.mAlphaDose{1}(V,:)*wOnes + (dij.mSqrtBetaDose{1}(V,:)*wOnes).^2);
           % ensure a underestimated biological effective dose 
           TolEstBio        = 1.2;
           % calculate maximal RBE in target
           maxCurrRBE = max(-cst{ixTarget,5}.alphaX + sqrt(cst{ixTarget,5}.alphaX^2 + ...
                        4*cst{ixTarget,5}.betaX.*CurrEffectTarget)./(2*cst{ixTarget,5}.betaX*(dij.physicalDose{1}(V,:)*wOnes)));
           wInit    =  ((doseTarget)/(TolEstBio*maxCurrRBE*max(dij.physicalDose{1}(V,:)*wOnes)))* wOnes;
    end
    
else 
    bixelWeight =  (doseTarget)/(mean(dij.physicalDose{1}(V,:)*wOnes)); 
    wInit       = wOnes * bixelWeight;
    pln.propOpt.bioOptimization = 'none';
end

% set optimization options
options.radMod          = pln.radiationMode;
options.bioOpt          = pln.propOpt.bioOptimization;
options.ID              = [pln.radiationMode '_' pln.propOpt.bioOptimization];
options.numOfScenarios  = dij.numOfScenarios;

%% Select Projections (i.e., quantities that are optimized)
switch pln.propOpt.bioOptimization
    case 'LEMIV_effect'
        backProjection = matRad_EffectProjection;
    case 'const_RBExD'
        backProjection = matRad_ConstantRBEProjection;
    case 'LEMIV_RBExD'
        backProjection = matRad_VariableRBEProjection;
    case 'none'
        backProjection = matRad_DoseProjection;
    otherwise
        warning(['Did not recognize bioloigcal setting ''' pln.probOpt.bioOptimization '''!\nUsing physical dose optimization!']);
        backProjection = matRad_DoseProjection;
end

optiProb = matRad_OptimizationProblem(backProjection);

%If LET objectives are present, also add an LET projection
if haveLETobjectives
    optiProb.BP_LET = matRad_LETProjection();
end

if haveDADRobjectives
    optiProb.BP_DADRfixed = matRad_DADRProjectionFixedCurrent();
end


%% Handle Optimizer
if ~isfield(pln.propOpt,'optimizer')
    pln.propOpt.optimizer = 'IPOPT';
end

switch pln.propOpt.optimizer
    case 'IPOPT'
        optimizer = matRad_OptimizerIPOPT;
    case 'fmincon'
        optimizer = matRad_OptimizerFmincon;
    otherwise
        warning(['Optimizer ''' pln.propOpt.optimizer ''' not known! Fallback to IPOPT!']);
        optimizer = matRad_OptimizerIPOPT;
end
        
%optimizer = matRad_OptimizerFmincon;

optimizer = optimizer.optimize(wInit,optiProb,dij,cst);

wOpt = optimizer.wResult;
info = optimizer.resultInfo;

resultGUI = matRad_calcCubes(wOpt,dij);
resultGUI.wUnsequenced = wOpt;
resultGUI.usedOptimizer = optimizer;
resultGUI.info = info;

% unblock mex files
clear mex
