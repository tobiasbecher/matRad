
matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 500000;
matRad_cfg.propOpt.defaultAccChangeTol = 1e-6;
%%
load('PROSTATE.mat');

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
pln(1).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.z = 5; % [mm]

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
pln(2).propStf.gantryAngles    = [0:52:359];%[0:90:359]; % [?] ;
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
pln(2).propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.z = 5; % [mm]
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
%%


pln_1 = matRad_plnWrapper(pln(1));

% Stf
stf_singleModalityOneJO = matRad_stfWrapper(ct,cst, pln_1);

% dij
dij_singleModalityOneJO = matRad_calcCombiDose(ct,stf_singleModalityOneJO,pln_1, cst, 0);


%%
pln_2 = matRad_plnWrapper(pln(2));

% Stf
stf_singleModalityTwoJO = matRad_stfWrapper(ct,cst, pln_2);

% dij
dij_singleModalityTwoJO = matRad_calcCombiDose(ct,stf_singleModalityTwoJO,pln_2, cst, 0);
%%
pln_M = matRad_plnWrapper(pln);

% Stf
stf_M = matRad_stfWrapper(ct,cst, pln_M);

% dij
%%     
dij_M = matRad_calcCombiDose(ct,stf_M, pln_M, cst, 0);

%%normalization for physical dose
factor = mean(dij_M.original_Dijs{2}.physicalDose{1},'all')/mean(dij_M.original_Dijs{1}.physicalDose{1},'all')
dij_M.original_Dijs{1}.physicalDose{1} = dij_M.original_Dijs{1}.physicalDose{1}*factor;
%%

%%
%% No fem heads



PriorityList1 = matRad_PriorityList1();
PriorityList1.addObjective(1,DoseObjectives.matRad_SquaredDeviation(100,68),2,6);
PriorityList1.addObjective(2,DoseObjectives.matRad_SquaredDeviation(100,56),2,7);
PriorityList1.addObjective(3,DoseObjectives.matRad_EUDMin(100,0,8),20,1);
PriorityList1.addObjective(4,DoseObjectives.matRad_EUDMin(100,0,3),30,8);
PriorityList1.addObjective(5,DoseObjectives.matRad_MeanDose(100,0),10,9);
%constraints
PriorityList1.addConstraint(DoseConstraints.matRad_MinMaxDose(0,58),9);
PriorityList1.addConstraint(DoseConstraints.matRad_MinMaxDose(0,58),1);
PriorityList1.addConstraint(DoseConstraints.matRad_MinMaxDose(0,58),8);

%%
cst{2,5}.Priority = 5;
cst{4,5}.Priority = 6;
cst{10,5}.Priority = 7;
cst{9,5}.Priority = 8;
%%
%{
tic;
[resultGUIs,resultGUIs2,cst1,cst2,PriorityList2] =  matRad_2pecOptimizationNew(PriorityList1,dij_singleModalityOneJO,ct,cst,pln_1); 
tProtons = toc;
%}
%%


%% Femoral heads with mixed Mods


PriorityListm = matRad_PriorityList1();
PriorityListm.addObjective(1,DoseObjectives.matRad_SquaredDeviation(100,68),2,6);
PriorityListm.addObjective(2,DoseObjectives.matRad_SquaredDeviation(100,56),2,7);
PriorityListm.addObjective(3,DoseObjectives.matRad_EUDMin(100,0,8),20,1);
PriorityListm.addObjective(4,DoseObjectives.matRad_EUDMin(100,0,3),30,8);
PriorityListm.addObjective(5,DoseObjectives.matRad_EUDMin(100,0),20,4);
PriorityListm.addObjective(5,DoseObjectives.matRad_EUDMin(100,0),20,10);
PriorityListm.addObjective(6,DoseObjectives.matRad_MeanDose(100,0),10,9);
%constraints
PriorityListm.addConstraint(DoseConstraints.matRad_MinMaxDose(0,58),9);
PriorityListm.addConstraint(DoseConstraints.matRad_MinMaxDose(0,58),1);
PriorityListm.addConstraint(DoseConstraints.matRad_MinMaxDose(0,58),8);
PriorityListm.addConstraint(DoseConstraints.matRad_MinMaxDose(0,58),4);
PriorityListm.addConstraint(DoseConstraints.matRad_MinMaxDose(0,58),10);
%%
cst{2,5}.Priority = 5;
cst{4,5}.Priority = 6;
cst{10,5}.Priority = 7;
cst{9,5}.Priority = 8;

%%
tic;
[resultGUIsM,resultGUIs2M,cst1M,cst2M,PriorityList2M] =  matRad_2pecOptimizationNew(PriorityListm,dij_M,ct,cst,pln_M); 
tMixed = toc;

%%


pln = pln(1)








%%
matRadGUI


























%%


%% calculate resultGUI scaled to dose grid
resultGUIProtonsDoseGrid = matRad_calcCubesDoseGridMixedWrapper(resultGUIs2{4},dij_singleModalityOneJO,1)
resultGUIProtons = matRad_calcCubesMixedWrapper(resultGUIs2{4},dij_singleModalityOneJO,1)
%% calculate indicator functions for mixed Plan on dose grid
cstR = cst;
cstR(10,:) = [];
cstR(5,:) = [];
cstR(4,:) = [];
cstR(3,:) = [];
cstR(2,:) = [];

%%
cstR = matRad_setOverlapPriorities(cstR)
cstR = matRad_resizeCstToGrid(cstR,dij_singleModalityOneJO.ctGrid.x,  dij_singleModalityOneJO.ctGrid.y,  dij_singleModalityOneJO.ctGrid.z,...
                                 dij_singleModalityOneJO.doseGrid.x,dij_singleModalityOneJO.doseGrid.y,dij_singleModalityOneJO.doseGrid.z);
matRad_indicatorWrapperMixed(matRad_setOverlapPriorities(cstR),pln_1,test,70)
%%

figure
dvh = matRad_calcDVH(cstR,test{1}.physicalDose*5)
matRad_showDVH(dvh,cstR,pln_1)

matlab2tikz('width','\fwidth','height','\fheight')
%%

%%
figure
matRad_plotSliceWrapper(gca(),ct,cstMod,1,test2{1}.physicalDose*5,3,slice)
zoom(1.4)

%%

matlab2tikz('width','\fwidth','height','\fheight')



%%
test = matRad_calcCubesDoseGridMixedWrapper(resultGUIs2M{6},dij_M,1)
test2 = matRad_calcCubesMixedWrapper(resultGUIs2M{6},dij_M,1)
%% calculate indicator functions for mixed Plan on dose grid
cstR = cst;
cstR(5,:) = [];
cstR(3,:) = [];
cstR(2,:) = [];

%%
cstR = matRad_setOverlapPriorities(cstR)
cstR = matRad_resizeCstToGrid(cstR,dij_singleModalityOneJO.ctGrid.x,  dij_singleModalityOneJO.ctGrid.y,  dij_singleModalityOneJO.ctGrid.z,...
                                 dij_singleModalityOneJO.doseGrid.x,dij_singleModalityOneJO.doseGrid.y,dij_singleModalityOneJO.doseGrid.z);
matRad_indicatorWrapperMixed(matRad_setOverlapPriorities(cstR),pln_M,test,70)
%%

figure
dvh = matRad_calcDVH(cstR,test{1}.physicalDose*5+test{2}.physicalDose*25)
matRad_showDVH(dvh,cstR,pln_1)
%%

matlab2tikz('width','\fwidth','height','\fheight')
%%
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
cstMod = cst;
cstMod(5,:) = [];
cstMod(3,:) = [];
cstMod(2,:) = [];
%%
figure
matRad_plotSliceWrapper(gca(),ct,cstMod,1,test2{1}.physicalDose*5+test2{2}.physicalDose*25,3,slice)
zoom(1.4)

%%
figure
matRad_plotSliceWrapper(gca(),ct,cstMod,1,test2{1}.physicalDose*5+test2{2}.physicalDose*25,3,slice)
zoom(1.4)

%%

matlab2tikz('width','\fwidth','height','\fheight')
%}
%%
matRad_plotLexicographicMixed(pln_M,ct,cst,resultGUIs2M)
%%
resultGUI = matRad_convertMixedPlan(resultGUIs2M{5},pln_M)
%%
pln = pln(1)
%%
matRadGUI


%%


























%% ANALYSIS PROTONS


resultGUIProtonsDoseFem = matRad_calcCubesDoseGridMixedWrapper(resultGUIs2{4},dij_singleModalityOneJO,1)
resultGUIProtonsDose = matRad_calcCubesMixedWrapper(resultGUIs2{4},dij_singleModalityOneJO,1)

cstR = cst;
cstR(5,:) = [];
cstR(3,:) = [];
cstR(2,:) = [];



%%
cstR = matRad_setOverlapPriorities(cstR)
cstR = matRad_resizeCstToGrid(cstR,dij_singleModalityOneJO.ctGrid.x,  dij_singleModalityOneJO.ctGrid.y,  dij_singleModalityOneJO.ctGrid.z,...
                                 dij_singleModalityOneJO.doseGrid.x,dij_singleModalityOneJO.doseGrid.y,dij_singleModalityOneJO.doseGrid.z);
matRad_indicatorWrapperMixed(matRad_setOverlapPriorities(cstR),pln_1,resultGUIProtonsDoseFem,70)

%%
figure
matRad_plotSliceWrapper(gca(),ct,cstR,1,resultGUIProtonsDose{1}.physicalDose*5,3,slice)
zoom(1.4)

%%
figure
dvhProtons = matRad_calcDVH(cstR,resultGUIProtonsDoseFem{1}.physicalDose*5)
matRad_showDVH(dvhProtons,cstR,pln_1)
dvhMixed = matRad_calcDVH(cstR,test{1}.physicalDose*5+test{2}.physicalDose*25)
matRad_showDVH(dvhMixed,cstR,pln_1,2)





