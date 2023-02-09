function returnStruct = matRad_paretoRennen(dij,cst,pln,VOIs,wInit)
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

%loop over VOIs
for  i = 1:numel(VOIStr)
    foundVOI = false;
    %loop over cst volumes
    for j = 1:size(cst,1)
        if VOIStr(i)== cst{j,2} %is it an objective we are interest in?
            foundVOI = true;
            %need to check if all are doseobjectives or if there are constraints
            idxVOI(i) = j;
            VOIObjCount = 0;
            %sizes(i) =  size(cst{j,6},2);
            
            for k = 1:size(cst{j,6},2)
                if  contains(class(cst{j,6}{k}),'DoseObjectives') % is it an objective or constraint?
                    name = VOIStr(i) + " " + convertCharsToStrings(cst{j,6}{k}.name);
                    VOIObjNames = [VOIObjNames name];
                    VOIObjCount = VOIObjCount+1;
                end
            end
            %sanitycheck to see if the VOI given actually contains objectives
            if VOIObjCount == 0
                
               matRad_cfg.dispError('Chosen VOI "%s" contains no Dose Objectives Please choose another one!\n',VOIStr(i));
            end
            sizes(i) = VOIObjCount;
        end 
    end
    if ~foundVOI
        matRad_cfg.dispError('Chosen VOI "%s" not found! (Check spelling and capital letters)\n',VOIStr(i));
    end
end


numOfObj = sum(sizes);
fprintf('NumOfObj: %d \n',numOfObj);

%% generaete Anchor Points for optimization
[pen,penGrid] = matRad_generateAnchorPoints(sizes);
matRad_plotPenaltyGrid(penGrid);

weights = zeros(numel(wInit),size(penGrid,1));
fInd = zeros(size(penGrid,1),numOfObj);

% loop over all penalty combinations
for i = 1:size(pen{1},1)
    fprintf('Iteration: %d \n',i);
    
    % loop over structures of interest and update penValues
    % !only loops over structures with varying penalties so far!
    for j = 1:numel(idxVOI) %loop over indices
        
        for k = 1:size(pen{j},2)
            if contains(class(cst{idxVOI(j),6}{k}),'DoseObjectives') % only consider objectives, not constraints
                cst{idxVOI(j),6}{k}.penalty = pen{j}(i,k)*100;
            end
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
    %update initial weights
    wInit = wOpt;
    %figure(fig1), plot(fInd{1},fInd{2});
    %refreshdata
    %drawnow
    pause(1)
end
time = toc;


%% Generaete further points

%initialize OPS boundaries
OPSA = [];
OPSb = [];

% first on: "Balanced point" after anchor points (slightly different calculation for normal)
%
np = size(penGrid,1);
[a,b,firstNormal] = matRad_normalFromFacet(penGrid,1:np,1);

penGrid(np+1,:) = abs(firstNormal);
    

%update cst values
pidx = 1;
for j = 1:numel(idxVOI) %loop over indices
    
    for k = 1:size(cst{idxVOI(j),6},2)
        if contains(class(cst{idxVOI(j),6}{k}),'DoseObjectives') % only consider objectives, not constraints
            cst{idxVOI(j),6}{k}.penalty = newPen(pidx)*100;
            pidx = pidx+1;
        end
    end
end

%calculate new point

optimizer = optimizer.optimize(wInit,optiProb,dij,cst);
wOpt = optimizer.wResult;
info = optimizer.resultInfo;
weights = [weights,wOpt];
wInit = wOpt;
fIndv = matRad_objectiveFunctions(optiProb,wOpt,dij,cst)
fInd = [fInd;fIndv]

OPSA = [OPSA;-penalties(np+1,:)];
OPSb = [OPSb;-1*ps(np+1,:)*penalties(np+1,:)'];

errors = [];






%% remaining facets
nIter = 50;
for i = 1:nIter
    %Step 1 calculate convex Hull -> Inner approximation (IPS) and gives facets
    %Rennen Algorithm
    fVals = fInd(1:size(penGrid,1)-nIter+i-1,:);
    
    %calculate epsilon value
    L = min(fVals,[],1);
    U = max(fVals,[],1);
    eps = U - L;


    fValsMod = matRad_generateDummyPoints(fVals); %generate dummy points
    %
    [k,vol] = convhulln(fValsMod);
    [kred,vol] = convhulln(fVals);
    %check for relevant facets (those that contain points of the original
    %fVals set)
    IPSidxs = 1:size(fVals,1);
    relFacetidxs = [];
            
    for j = 1:size(k,1)
        if any(ismember(k(j,:),IPSidxs))
            relFacetidxs = [relFacetidxs,j];
        end
    end
    facetMods = k(relFacetidxs,:);
    facetErrors = zeros(size(facetMods,1),1);
    normals = zeros(size(facetMods));
    
    %% 

    %Loop over facets ands 
    for j = 1:size(relFacetidxs,2)
        [facetPoints,refPoint,normal] = matRad_normalFromFacet(fValsMod,facetMods,j);

        
        %check for sign of normals (should be redundant)
        if all(normal<0)
            continue
        end
        
        %now check for OPS point for facet
        lb = min(fVals,[],1);
        ub = max(fVals,[],1);
        z = linprog(normal,OPSA,OPSb,[],[],lb,ub); 
        
        %hyperplane distance
        b = refPoint*normal;

        %calculate error for each facet
        
        facetErrors(j) = (b-z'*normal)/(eps*normal); 
        normals(j,:) = normal;
        %{
        figure
        trisurf(k,fValsMod(:,1),fValsMod(:,2),fValsMod(:,3),'FaceColor','cyan')
        hold on 
        fill3(facetPoints(:,1),facetPoints(:,2),facetPoints(:,3),'green')
        %}
    end

    [A,I] = sort(facetErrors,'descend');

    %%check for next facet to run
    found = false;
    facetNum= 1;    
    w = zeros(1,size(penalties,2));
    accuracy = 3;

    while ~found && facetNum <= numel(I) %loop over facets
        idx = I(facetNum);
        norm = normals(idx,:);
        if ~any(ismember(round(penalties,accuracy),round(norm,accuracy),'rows'))
            errors = [errors,facetErrors(idx)];
            newPen = norm;
            %update weights in cst
            pidx = 1;
            for j = 1:numel(idxVOI) %loop over indices
                
                for k = 1:size(cst{idxVOI(j),6},2)
                    if contains(class(cst{idxVOI(j),6}{k}),'DoseObjectives') % only consider objectives, not constraints
                        cst{idxVOI(j),6}{k}.penalty = newPen(pidx)*100;
                        pidx = pidx+1;
                    end
                end
            end
            
            optimizer = optimizer.optimize(wInit,optiProb,dij,cst);
            wOpt = optimizer.wResult;
            info = optimizer.resultInfo;
            weights = [weights,wOpt];
            penGrid = [penGrid;newPen];
            
            fIndv = matRad_objectiveFunctions(optiProb,wOpt,dij,cst)
            fInd = [fInd;fIndv]
            fprintf('fInd%d\n',fInd);


            found = true;
        end
        facetNum = facetNum +1;
    end    


    % when final point is found: Update OPsw and OPSb
    OPSA = [OPSA;-w]; %add normal vector of facet that was run 
    OPSb = [OPSb;-ps(i+4,:)*w'];

end





returnStruct.VOIObj = VOIObjNames;
returnStruct.pen = pen;
returnStruct.time = time;
returnStruct.weights = weights;
returnStruct.finds = fInd;
returnStruct.finds2 = fInd2;
returnStruct.penGrid = penGrid;
% unblock mex files
clear mex
