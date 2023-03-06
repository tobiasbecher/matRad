function returnStruct = matRad_RennenRadio(dij,cst,pln,wInit)
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
%get number of objectives at the start
objcount = 0;
for i = 1:size(cst,1) % loop over cst
    for j = 1:numel(cst{i,6})
        %check whether dose objective or constraint (should all be class
        %here no structs)
        if contains(class(cst{i,6}{j}),'DoseObjectives')
            objcount = objcount + 1; %update current objective number
        end
    end
end

%% generaete Anchor Points for optimization
[pen,penGrid] = matRad_generateAnchorPoints(ones(objcount,1));
matRad_plotPenaltyGrid(penGrid);

weights = zeros(numel(wInit),size(penGrid,1));
fInd = zeros(size(penGrid,1),objcount);
objectiveFunctionVals = {};
% loop over all penalty combinations

for i = 1:size(pen{1},1)
    cst = matRad_updatecst(cst,penGrid(i,:));
    
    
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
    
    fInd(i,:) = matRad_objectiveFunctions(optiProb,wOpt,dij,cst);
    objectiveFunctionVals(end + 1) = {optimizer.allObjectiveFunctionValues};
    
end



%% try to normalize objectives

optiProb.normalizationScheme.type = 'UL'
optiProb.normalizationScheme.U = max(fInd,[],1);
optiProb.normalizationScheme.L = min(fInd,[],1);

fInd = (fInd-optiProb.normalizationScheme.L)./(optiProb.normalizationScheme.U-optiProb.normalizationScheme.L)

'HI'
% loop over all penalty combinations






%{
varToLoad = {'data'}
load('dataRennen.mat',varToLoad{:})
weights = data.weights;
penGrid = data.penGrid(1:4,:);
weights = weights(:,1:4);
fInd = data.finds;
fInd = fInd(1:4,:);
%}

%% Grouping
if isfield(pln.propOpt,'grouping') %if no grouping ignore
    
    switch pln.propOpt.grouping
        case 'forced' %grouping of objectives forced by user -> How to store groups? Pln para?
            if isfield(pln.propOpt,'groups')
                groups = pln.propOpt.groups; %?  
            else
                 matRad_cfg.dispError('pln.propOpt.grouping set to forced, but no groups provided!');
            end
        case 'Spearman' %grouping by correlation of spearman ranking
                %tbd
            'Spearman ranking: TO BE IMPLEMENTED'

    end
end

%Do we have to reoptimize as fInd would change size? 
%{
objcount = numel(groups); %update 

fInd = zeros(objcount);
[pen,penGrid] = matRad_generateAnchorPoints(objcount);
weights = zeros(numel(wInit),objcount);
%}
%% Generaete further points

%initialize OPS boundaries
OPSA = [];
OPSb = [];
np = size(penGrid,1);

% first on: "Balanced point" after anchor points (slightly different calculation for normal)
%

[a,b,firstNormal] = matRad_normalFromFacet(fInd,1:np,1);

penGrid(np+1,:) = abs(firstNormal);
newPen = abs(firstNormal);

%update cst values

cst = matRad_updatecst(cst,newPen); % change for groupings!

%calculate new point

optimizer = optimizer.optimize(wInit,optiProb,dij,cst);
wOpt = optimizer.wResult;
info = optimizer.resultInfo;
weights = [weights,wOpt];
objectiveFunctionVals(end + 1) = {optimizer.allObjectiveFunctionValues};
%wInit = wOpt;
fIndv = matRad_objectiveFunctions(optiProb,wOpt,dij,cst);
fIndv = optiProb.normalizeObjectives(fIndv);
%fIndv = cellfun(@(group)sum(fIndv(:,group)),groups)
%grouping

%change for groupings -> sum up over groups!
fInd = [fInd;fIndv];

removedfInd = [];
removedpenGrid = [];
removedweights = [];
removedOPSA = [];
removedOPSb = [];


OPSA = [OPSA;-penGrid(np+1,:)];
OPSb = [OPSb;-1*fInd(np+1,:)*penGrid(np+1,:)'];
errors = [];
allErrors = {};

initSize = size(penGrid,1);




%% remaining facets
nIter = 50;
for i = 1:nIter
    fprintf('Now in iteration %i',i)
    %Step 1 calculate convex Hull -> Inner approximation (IPS) and gives facets
    %Rennen Algorithm
    fVals = fInd;
    
    %calculate epsilon value (used in error calculation)
    L = min(fVals,[],1);
    U = max(fVals,[],1);
    eps = U - L;


    fValsMod = matRad_generateDummyPoints(fVals); %generate dummy points
    %
    [kmod,vol] = convhulln(fValsMod);
    [kred,vol] = convhulln(fVals);
    %check for relevant facets (those that contain points of the original
    %fVals set)
    IPSidxs = 1:size(fVals,1);
    relFacetidxs = [];
            
    for j = 1:size(kmod,1)
        if any(ismember(kmod(j,:),IPSidxs))
            relFacetidxs = [relFacetidxs,j];
        end
    end
    facetMods = kmod(relFacetidxs,:);
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
        trisurf(kmod,fValsMod(:,1),fValsMod(:,2),fValsMod(:,3),'FaceColor','cyan')
        hold on 
        fill3(facetPoints(:,1),facetPoints(:,2),facetPoints(:,3),'green')
        %}
    end
    allErrors(end+1) = {facetErrors};
    [A,I] = sort(facetErrors,'descend');

    %%check for next facet to run
    found = false;
    facetNum= 1;
    accuracy = 4;

    while ~found && facetNum <= numel(I) %loop over facets
        idx = I(facetNum);
        norm = normals(idx,:);
        %check if facet has already been run
        if ~any(ismember(round(penGrid,accuracy),round(norm,accuracy),'rows'))

            %update weights in cst
            
            cst = matRad_updatecst(cst,norm);
            paretoOptimal = false;
            calcPointsForFacet = 0;
            while ~paretoOptimal && calcPointsForFacet < 4 %more than 4 points?
                if calcPointsForFacet < 3
                    optimizer = optimizer.optimize(weights(:,end-calcPointsForFacet),optiProb,dij,cst);
                else
                    optimizer = optimizer.optimize(wInit,optiProb,dij,cst);
                end
                %really necessary?
                if optimizer.resultInfo.iter < 15 %sometimes optimizer will terminate fast (points tend to be not optimal!)
                    calcPointsForFacet = calcPointsForFacet + 1;
                    continue
                end
                %check number of iterations needed for optimzation
                wOpt = optimizer.wResult;
                
                fIndv = matRad_objectiveFunctions(optiProb,wOpt,dij,cst);
                fIndv = optiProb.normalizeObjectives(fIndv);
                %how does the newly generated point influence the pareto
                %surface?
    
                dominatedPoints = matRad_paretoCheck(fInd,fIndv); 
                % if no point is dominated -> returns []
                % if the newly generated point is not optimal -> Returns array with a 0 inside 
                % if (a) previously generated point(s) is/are dominated by the
                % newly generated point -> returns index in fInd -> Need to
                % update penGrid, weights, fInd,OPSa,OPSb

                if isempty(dominatedPoints)
                    
                    paretoOptimal = true;   
                    
                elseif dominatedPoints == 0
                    calcPointsForFacet = calcPointsForFacet + 1;
                    continue %% found point is dominated by previously calculated point -> has to be recalculated

                else %newly calculated point dominates points on the previous pareto surface
                    paretoOptimal = true;
                    %remove dominated points and related equations %might
                    %have to check if this fully works
                    removedfInd  = [removedfInd;fInd(dominatedPoints,:)];
                    removedpenGrid = [removedpenGrid;penGrid(dominatedPoints,:)];
                    removedweights = [removedweights,weights(:,dominatedPoints)];


                    fInd(dominatedPoints,:) = []; 
                    penGrid(dominatedPoints,:) = []; 
                    weights(:,dominatedPoints) = [];
                    %Issue if 
                    if (dominatedPoints-objcount <=0)
                      'Well well well'
                    else
                        removedOPSA = [removedOPSA;OPSA(dominatedPoints-objcount,:)];
                        removedOPSb = [removedOPSb;OPSb(dominatedPoints-objcount,:)];
                        OPSA(dominatedPoints-objcount,:) = []; % shifted index
                        OPSb(dominatedPoints-objcount,:) = []; 

                    end
                    %remove dominated points
                    
                end
            end
            objectiveFunctionVals(end + 1) = {optimizer.allObjectiveFunctionValues};
            errors = [errors,facetErrors(idx)];
            newPen = norm;
            %awInit = wOpt;
            info = optimizer.resultInfo;
            weights = [weights,wOpt];
            penGrid = [penGrid;newPen];

            
            %fIndv = cellfun(@(group)sum(fIndv(:,group)),groups)
            fInd = [fInd;fIndv];
            
            %%  plot current convex hull and facet that was run            
               
            figure
            trisurf(kmod,fValsMod(:,1),fValsMod(:,2),fValsMod(:,3),'FaceColor','cyan')
            hold on
            fill3(fValsMod(facetMods(idx,:),1),fValsMod(facetMods(idx,:),2),fValsMod(facetMods(idx,:),3),'r')

            found = true;
        end
        facetNum = facetNum +1;
    end    


    % when final point is found: Update OPSA and OPSb

    OPSA = [OPSA;-newPen]; %add normal vector of facet that was run 
    OPSb = [OPSb;-fIndv*newPen'];

end





%returnStruct.VOIObj = VOIObjNames;
returnStruct.OPSA = OPSA;
returnStruct.OPSb = OPSb;
returnStruct.pen = pen;
returnStruct.weights = weights;
returnStruct.finds = fInd;
returnStruct.penGrid = penGrid;
returnStruct.errors = errors;
returnStruct.allErrors = allErrors;
returnStruct.removed = {removedfInd,removedpenGrid,removedweights,removedOPSA,removedOPSb};
% unblock mex files
%}
returnStruct.weights = weights;
returnStruct.finds = fInd;
returnStruct.penGrid = penGrid;


returnStruct.allObj = objectiveFunctionVals;
clear mex
