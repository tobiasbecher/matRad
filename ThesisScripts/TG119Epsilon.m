
matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 500000;
matRad_cfg.propOpt.defaultAccChangeTol = 1e-6/30;
%%
load('TG119.mat');

%%

pln(1).numOfFractions  = 5;
pln(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln(1).machine         = 'Generic';

% beam geometry settings
pln(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(1).propStf.gantryAngles    = [90 270];
pln(1).propStf.couchAngles     = zeros(numel(pln(1).propStf.gantryAngles),1); % [?] ; 
pln(1).propStf.numOfBeams      = numel(pln(1).propStf.gantryAngles);
pln(1).propStf.isoCenter       = ones(pln(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln(1).propDoseCalc.calcLET = 1;

pln(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(1).propOpt.spatioTemp      = 0;
pln(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln(1).propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.z = 2.5; % [mm]

% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD

%=======================================> Model check error in bioModel
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln(1).bioParam = matRad_bioModel(pln(1).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln(1).multScen = matRad_multScen(ct,scenGenType);

%% Define meta info for pln(2) 
pln(2).numOfFractions  = 25;
pln(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln(2).machine         = 'Generic';

% beam geometry settings
pln(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(2).propStf.gantryAngles    = [0:40:359];%[0:90:359]; % [?] ;
pln(2).propStf.couchAngles     = zeros(numel(pln(2).propStf.gantryAngles),1);  % [?] ; 
pln(2).propStf.numOfBeams      = numel(pln(2).propStf.gantryAngles);
pln(2).propStf.isoCenter       = ones(pln(2).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln(2).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(2).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(2).propOpt.spatioTemp      = 0;
pln(2).propOpt.STscenarios     = 5;
%pln(2).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln(2).propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.z = 2.5; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln(2).bioParam = matRad_bioModel(pln(2).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln(2).multScen = matRad_multScen(ct,scenGenType);

% prepping cst 
% placing alpha/beta ratios in cst{:,6},
% different alpha beta ration for each obj of a structure  
sparecst = 0;
%%
cst = matRad_prepCst(cst, sparecst);

%% Plan definitions

pln_2 = matRad_plnWrapper(pln(2));

% Stf
stf_singleModalityTwoJO = matRad_stfWrapper(ct,cst, pln_2);

% dij
dij_singleModalityTwoJO = matRad_calcCombiDose(ct,stf_singleModalityTwoJO,pln_2, cst, 0);



%% Definition of priority list
PriorityList1 = matRad_PriorityList1();
PriorityList1.addObjective(1,DoseObjectives.matRad_SquaredDeviation(100,50),1.5,2);
PriorityList1.addObjective(2,DoseObjectives.matRad_EUDMin(100,0),22.5,1);
PriorityList1.addObjective(3,DoseObjectives.matRad_MeanDose(100,0),7,3);
%constraints

PriorityList1.addConstraint(DoseConstraints.matRad_MinMaxDose(0,40),1);
PriorityList1.addConstraint(DoseConstraints.matRad_MinMaxDose(45,57),2);
PriorityList1.addConstraint(DoseConstraints.matRad_MinMaxDose(0,45),3);

%% Calculate plan using 2p method

tic;
[resultGUIs,resultGUIs2,cst1,cst2,PriorityList2] =  matRad_2pecOptimizationNew(PriorityList1,dij_singleModalityTwoJO,ct,cst,pln_2); 
tNoFemHeads = toc;
%%
PriorityList1 = matRad_PriorityList1();
PriorityList1.addObjective(1,DoseObjectives.matRad_SquaredDeviation(100,50),1.3,2);
PriorityList1.addObjective(2,DoseObjectives.matRad_EUDMin(100,0),15,1);
PriorityList1.addObjective(3,DoseObjectives.matRad_MeanDose(100,0),7,3);
%constraints

PriorityList1.addConstraint(DoseConstraints.matRad_MinMaxDose(0,40),1);
PriorityList1.addConstraint(DoseConstraints.matRad_MinMaxDose(45,57),2);
PriorityList1.addConstraint(DoseConstraints.matRad_MinMaxDose(0,45),3);

%% Calculate plan using 2p method

tic;
[resultGUIs,resultGUIs2,cst1,cst2,PriorityList2] =  matRad_2pecOptimizationNew(PriorityList1,dij_singleModalityTwoJO,ct,cst,pln_2); 
tNoFemHeads = toc;



%%
PriorityList1 = matRad_PriorityList1();
PriorityList1.addObjective(1,DoseObjectives.matRad_SquaredDeviation(100,50),2,2);
PriorityList1.addObjective(2,DoseObjectives.matRad_EUDMin(100,0),10,1);
PriorityList1.addObjective(3,DoseObjectives.matRad_MeanDose(100,0),7,3);
%constraints

PriorityList1.addConstraint(DoseConstraints.matRad_MinMaxDose(0,40),1);
PriorityList1.addConstraint(DoseConstraints.matRad_MinMaxDose(45,57),2);
PriorityList1.addConstraint(DoseConstraints.matRad_MinMaxDose(0,45),3);

%% Calculate plan using 2p method
tic;
[resultGUIs,resultGUIs2,cst1,cst2,PriorityList2] =  matRad_2pecOptimizationNew(PriorityList1,dij_singleModalityTwoJO,ct,cst,pln_2); 
tNoFemHeads = toc;





%%
%load('DataCorrectObjAll3ConstrHighRes.mat')

%%
%resultGUI = matRad_convertMixedPlan(resultGUIs2{3},pln_2,82)
%%
resultGUI = matRad_convert2pPlans(resultGUIs,pln_2)

%%
matRadGUI

%% Analysis of plans


%calculate cubes
resultGUIPhotonsDoseGrid1 = matRad_calcCubesDoseGridMixedWrapper(resultGUIs2{1},dij_singleModalityTwoJO,1);
resultGUIPhotons1 = matRad_calcCubesMixedWrapper(resultGUIs2{1},dij_singleModalityTwoJO,1);
%%
resultGUIPhotonsDoseGrid2 = matRad_calcCubesDoseGridMixedWrapper(resultGUIs2{2},dij_singleModalityTwoJO,1);
resultGUIPhotons2 = matRad_calcCubesMixedWrapper(resultGUIs2{2},dij_singleModalityTwoJO,1);
%%

cstR = matRad_setOverlapPriorities(cst)
cstR = matRad_resizeCstToGrid(cstR,dij_singleModalityTwoJO.ctGrid.x,  dij_singleModalityTwoJO.ctGrid.y,  dij_singleModalityTwoJO.ctGrid.z,...
                                 dij_singleModalityTwoJO.doseGrid.x,dij_singleModalityTwoJO.doseGrid.y,dij_singleModalityTwoJO.doseGrid.z);
matRad_indicatorWrapperMixed(matRad_setOverlapPriorities(cstR),pln_2,resultGUIPhotonsDoseGrid1)
%% Calculate dvh


figure
dvh1 = matRad_calcDVH(cstR,resultGUIPhotonsDoseGrid1{1}.physicalDose*5)
dvh2 = matRad_calcDVH(cstR,resultGUIPhotonsDoseGrid2{1}.physicalDose*5)
matRad_showDVH(dvh1,cstR,pln_2)
matRad_showDVH(dvh2,cstR,pln_2,2)
%% 
matRad_plotLexicographicMixed(pln_2,ct,cst,resultGUIs2Point2)
%%
%matlab2tikz('width','\fwidth','height','\fheight')
%% Plot a slice of the final plan

slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);

%%

[hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUIPhotons2{1}.physicalDose*25,3,slice)
for j= 1: numel(hContour)
    hContour{j}(1).LineWidth = 1;
end
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'ylabel',[])
set(gca,'xlabel',[])
set(gca,'title',[])
zoom(1.4)
%%

matlab2tikz('width','\fwidth','height','\fheight')





