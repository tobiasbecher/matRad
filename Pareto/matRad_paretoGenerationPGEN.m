function returnStruct = matRad_paretoGenerationPGEN(dij,cst,pln,VOIs,wInit)
% matRad inverse pareto planning wrapper function
% 
% call
%   [resultGUI,optimizer] = matRad_paretoGeneration(dij,cst,pln,pen)
%   [resultGUI,optimizer] = matRad_paretoGeneration(dij,cst,pln,pen,wInit)
%
% input
%   dij:        matRad dij struct
%   cst:        matRad cst struct
%   pln:        matRad pln struct
%   nPoints:    Number of pareto optimal points
%   VOIs:       Volumes for variation of penalties
%   wInit:      (optional) custom weights to initialize problems
%
% output
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%   optimizer:  Used Optimizer Object
%   fInd:       Array storing each individual final obj fnct value
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

% consider VOI priorities
cst  = matRad_setOverlapPriorities(cst);

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
        if ~isa(obj,'matRad_DoseOptimizationFunction')
            try
                obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
            catch
                matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!',i,j);
            end
        end
        
        obj = obj.setDoseParameters(obj.getDoseParameters()/pln.numOfFractions);
        
        cst{i,6}{j} = obj;        
    end
end

% resizing cst to dose cube resolution 
cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,  dij.ctGrid.y,  dij.ctGrid.z,...
                                 dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);

% find target indices and described dose(s) for weight vector
% initialization
V          = [];
doseTarget = [];
ixTarget   = [];

for i = 1:size(cst,1)
    if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
        V = [V;cst{i,4}{1}];
        
        %Iterate through objectives/constraints
        fDoses = [];
        for fObjCell = cst{i,6}
            dParams = fObjCell{1}.getDoseParameters();
            %Don't care for Inf constraints
            dParams = dParams(isfinite(dParams));
            %Add to dose list
            fDoses = [fDoses dParams];
        end
                
        doseTarget = [doseTarget fDoses];
        ixTarget   = [ixTarget i*ones(1,length(fDoses))];
    end
end
[doseTarget,i] = max(doseTarget);
ixTarget       = ixTarget(i);
wOnes          = ones(dij.totalNumOfBixels,1);
   
% calculate initial beam intensities wInit
matRad_cfg.dispInfo('Estimating initial weights... ');
if exist('wInit','var')
    %do nothing as wInit was passed to the function
    matRad_cfg.dispInfo('chosen provided wInit!\n');   
elseif strcmp(pln.bioParam.model,'constRBE') && strcmp(pln.radiationMode,'protons')
    % check if a constant RBE is defined - if not use 1.1
    if ~isfield(dij,'RBE')
        dij.RBE = 1.1;
    end

    doseTmp = dij.physicalDose{1}*wOnes;
    bixelWeight =  (doseTarget)/(dij.RBE * mean(doseTmp(V)));     
    wInit       = wOnes * bixelWeight;
    matRad_cfg.dispInfo('chosen uniform weight of %f!\n',bixelWeight);  
        
elseif pln.bioParam.bioOpt
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
                matRad_cfg.dispWarning('Reference dose > 10 Gy[RBE] for target. Biological optimization outside the valid domain of the base data. Reduce dose prescription or use more fractions.\n');
            end
            
        end
    end
    
    dij.ixDose  = dij.bx~=0; 
        
    if isequal(pln.bioParam.quantityOpt,'effect')

           effectTarget = cst{ixTarget,5}.alphaX * doseTarget + cst{ixTarget,5}.betaX * doseTarget^2;
           aTmp = dij.mAlphaDose{1}*wOnes;
           bTmp = dij.mSqrtBetaDose{1} * wOnes;
           p = sum(aTmp(V)) / sum(bTmp(V).^2);
           q = -(effectTarget * length(V)) / sum(bTmp(V).^2);
           
           wInit        = -(p/2) + sqrt((p^2)/4 -q) * wOnes;

    elseif isequal(pln.bioParam.quantityOpt,'RBExD')

           %pre-calculations
           dij.gamma             = zeros(dij.doseGrid.numOfVoxels,1);   
           dij.gamma(dij.ixDose) = dij.ax(dij.ixDose)./(2*dij.bx(dij.ixDose)); 

            
           % calculate current effect in target
           aTmp = dij.mAlphaDose{1}*wOnes;
           bTmp = dij.mSqrtBetaDose{1} * wOnes;
           doseTmp = dij.physicalDose{1}*wOnes;

           CurrEffectTarget = aTmp(V) + bTmp(V).^2;
           % ensure a underestimated biological effective dose 
           TolEstBio        = 1.2;
           % calculate maximal RBE in target
           maxCurrRBE = max(-cst{ixTarget,5}.alphaX + sqrt(cst{ixTarget,5}.alphaX^2 + ...
                        4*cst{ixTarget,5}.betaX.*CurrEffectTarget)./(2*cst{ixTarget,5}.betaX*doseTmp(V)));
           wInit    =  ((doseTarget)/(TolEstBio*maxCurrRBE*max(doseTmp(V))))* wOnes;
    end
    matRad_cfg.dispInfo('chosen weights adapted to biological dose calculation!\n'); 
else 
    doseTmp = dij.physicalDose{1}*wOnes;
    bixelWeight =  (doseTarget)/mean(doseTmp(V));
    wInit       = wOnes * bixelWeight;
    matRad_cfg.dispInfo('chosen uniform weight of %f!\n',bixelWeight);   
end
    


%% calculate probabilistic quantities for probabilistic optimization if at least
% one robust objective is defined 

%Check how to use 4D data
if isfield(pln,'propOpt') && isfield(pln.propOpt,'scen4D')
    scen4D = pln.propOpt.scen4D;
else
    scen4D = 1; %Use only first 4D scenario for optimization
end

%If "all" provided, use all scenarios
if isequal(scen4D,'all')
    scen4D = 1:size(dij.physicalDose,1);
end

linIxDIJ = find(~cellfun(@isempty,dij.physicalDose(scen4D,:,:)))';

FLAG_CALC_PROB = false;
FLAG_ROB_OPT   = false;


for i = 1:size(cst,1)
    for j = 1:numel(cst{i,6})
       if strcmp(cst{i,6}{j}.robustness,'PROB') && numel(linIxDIJ) > 1
          FLAG_CALC_PROB = true;
       end
       if ~strcmp(cst{i,6}{j}.robustness,'none') && numel(linIxDIJ) > 1
          FLAG_ROB_OPT = true;
       end
    end
end

if FLAG_CALC_PROB
    [dij] = matRad_calculateProbabilisticQuantities(dij,cst,pln);
end


% set optimization options
if ~FLAG_ROB_OPT || FLAG_CALC_PROB     % if multiple robust objectives are defined for one structure then remove FLAG_CALC_PROB from the if clause
   ixForOpt = scen4D;
else
   ixForOpt = linIxDIJ;
end

switch pln.bioParam.quantityOpt
    case 'effect'
        backProjection = matRad_EffectProjection;
    case 'RBExD'
        %Capture special case of constant RBE
        if strcmp(pln.bioParam.model,'constRBE')
            backProjection = matRad_ConstantRBEProjection;
        else
            backProjection = matRad_VariableRBEProjection;
        end
    case 'physicalDose'
        backProjection = matRad_DoseProjection;
    otherwise
        warning(['Did not recognize bioloigcal setting ''' pln.probOpt.bioOptimization '''!\nUsing physical dose optimization!']);
        backProjection = matRad_DoseProjection;
end

%Give scenarios used for optimization
backProjection.scenarios    = ixForOpt;
backProjection.scenarioProb = pln.multScen.scenProb;

optiProb = matRad_OptimizationProblem(backProjection);
optiProb.quantityOpt = pln.bioParam.quantityOpt;
if isfield(pln,'propOpt') && isfield(pln.propOpt,'useLogSumExpForRobOpt')
    optiProb.useLogSumExpForRobOpt = pln.propOpt.useLogSumExpForRobOpt;
end

%Get Bounds
if ~isfield(pln.propOpt,'boundMU')
    pln.propOpt.boundMU = false;
end 

if pln.propOpt.boundMU
    if (isfield(dij,'minMU') || isfield(dij,'maxMU')) && ~isfield(dij,'numParticlesPerMU')
        matRad_cfg.dispWarning('Requested MU bounds but number of particles per MU not set! Bounds will not be enforced and standard [0,Inf] will be used instead!');
    elseif ~isfield(dij,'minMU') && ~isfield(dij,'maxMU')
        matRad_cfg.dispWarning('Requested MU bounds but machine bounds not defined in dij.minMU & dij.maxMU! Bounds will not be enforced and standard [0,Inf] will be used instead!');
    else
        if isfield(dij,'minMU')
            optiProb.minimumW = dij.numParticlesPerMU .* dij.minMU / 1e6;
            matRad_cfg.dispInfo('Using lower MU bounds provided in dij!\n')
        end

        if isfield(dij,'maxMU')
            optiProb.maximumW = dij.numParticlesPerMU .* dij.maxMU / 1e6;
            matRad_cfg.dispInfo('Using upper MU bounds provided in dij!\n')
        end
    end
else
    matRad_cfg.dispInfo('Using standard MU bounds of [0,Inf]!\n')
end

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

%check that number of penalties equals number of constraints


% PARETO PART
%loop over VOI and get indices in cst file
tic
sizes = zeros(1,numel(VOIs));
idxVOI = zeros(size(VOIs));
VOIStr = convertCharsToStrings(VOIs);
VOIObjNames = [];

%loop over VOI and associate to index in cst. Also Get number of obj
%functions 
for  i = 1:numel(VOIStr)
    for j = 1:size(cst,1)
        if VOIStr(i)== cst{j,2}
            idxVOI(i) = j;
            sizes(i) =  size(cst{j,6},2);
            
            for k = 1:size(cst{j,6},2)
                name = VOIStr(i) + " " + convertCharsToStrings(cst{j,6}{k}.name);
                VOIObjNames = [VOIObjNames name];
            end
        end 
    end
end

numOfObj = sum(sizes);
fprintf('NumOfObj: %d \n',numOfObj);
[pen,penGrid] = matRad_generatePGENGrid(sizes);
matRad_plotPenaltyGrid(penGrid);

weights = zeros(numel(wInit),size(penGrid,1));
fInd = zeros(size(penGrid,1),numOfObj);
fInd2 = zeros(size(penGrid,1),numOfObj);
% loop over all penalty combinations
for i = 1:size(pen{1},1)
    fprintf('Iteration: %d \n',i);
    
    % loop over structures of interest and update penValues
    % !only loops over structures with varying penalties so far!
    for j = 1:numel(idxVOI) %loop over indices
        
        for k = 1:size(pen{j},2)
            cst{idxVOI(j),6}{k}.penalty;
            cst{idxVOI(j),6}{k}.penalty = pen{j}(i,k)*100;
        end
        
    end
    
    
    optimizer = optimizer.optimize(wInit,optiProb,dij,cst);
    wOpt = optimizer.wResult;
    info = optimizer.resultInfo;
    weights(:,i) = wOpt;
    info
    %set values for warm start
    %optimizer.optionsWarmStart.use      = true;
    %optimizer.optionsWarmStart.zl       = info.zl;
    %optimizer.optionsWarmStart.zu       = info.zu;
    %optimizer.optionsWarmStart.lambda   = info.lambda;
    

    cst = matRad_individualObjectiveFunction(optiProb,wOpt,dij,cst);
    fID  = 1;
    %get the indivdual objective function values and return them all
    for j = 1:numel(idxVOI)
        for k = 1:size(pen{j},2)
            fInd(i,fID) = cst{idxVOI(j),6}{k}.objValue;
 
            fID = fID + 1;
        end
    end
    fIndv = matRad_objectiveFunctions(optiProb,wOpt,dij,cst)
    fInd2(i,:) = fIndv
    %update initial weights
    wInit = wOpt;
    %figure(fig1), plot(fInd{1},fInd{2});
    %refreshdata
    %drawnow
    pause(1)
end
time = toc;
%returnStruct.weights = weights;
%returnStruct.finds = fInd;
%returnStruct.finds2 = fInd2;
returnStruct.VOIObj = VOIObjNames;
returnStruct.penGrid = penGrid;
returnStruct.pen = pen;
returnStruct.time = time;
%% Generaete further points
figure
matRad_plotParetoSurface(fInd,penGrid,VOIObjNames);
for i=1:10 %temporary for loop TODO: Use while loop that terminates after hitting distance (tbd)
 
    [a,b,c,d,dists,newPen] = matRad_convexHull(fInd,penGrid);
    newPen
    penGrid = [penGrid;newPen];
    matRad_plotPenaltyGrid(penGrid);
    
    %update penalties TODO: NOT CONFORM WITH CELL ARRAY STRUCTURE OF
    %PENGRID
    pidx = 1;
    for j = 1:numel(idxVOI) %loop over indices
        
        for k = 1:size(cst{idxVOI(j),6},2)
            ff = newPen(pidx);
            ff
            cst{idxVOI(j),6}{k}.penalty = newPen(pidx)*100;
            hh = cst{idxVOI(j),6}{k}.penalty;
            pidx = pidx+1;
        end
    end
    
    optimizer = optimizer.optimize(wInit,optiProb,dij,cst);
    wOpt = optimizer.wResult;
    info = optimizer.resultInfo;
    weights = [weights,wOpt];
    wInit = wOpt;
    fIndv = matRad_objectiveFunctions(optiProb,wOpt,dij,cst)
    fInd = [fInd;fIndv]
    fprintf('fInd%d\n',fInd);
    
    %warm start block
    optimizer.optionsWarmStart.use      = true;
    optimizer.optionsWarmStart.zl       = info.zl;
    optimizer.optionsWarmStart.zu       = info.zu;
    optimizer.optionsWarmStart.lambda   = info.lambda;
    %matRad_plotParetoSurface(fInd,penGrid,VOIObjNames);
end

returnStruct.weights = weights;
returnStruct.finds = fInd;
returnStruct.finds2 = fInd2;
returnStruct.penGrid = penGrid;
% unblock mex files
clear mex
