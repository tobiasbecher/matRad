%% Load patient data
matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 50000;
matRad_cfg.propOpt.defaultAccChangeTol = 1e-2;
load('PROSTATE.mat')
%%
%% Define meta info for pln(1)
%matRadGUI
%%
pln(1).numOfFractions  = 5;
pln(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln(1).machine         = 'Generic';

% beam geometry settings
pln(1).propStf.bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
pln(1).propStf.gantryAngles    = [0,180];
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
pln(2).numOfFractions  = 15;
pln(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln(2).machine         = 'Generic';

% beam geometry settings
pln(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(2).propStf.gantryAngles    = [0:50:359];%[0:90:359]; % [?] ;
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

%% Prepare cst

cst = matRad_prepCst(cst, sparecst);


%%

pln_2 = matRad_plnWrapper(pln(2));

% Stf
stf_singleModalityTwoJO = matRad_stfWrapper(ct,cst, pln_2);

% dij
dij_singleModalityTwoJO = matRad_calcCombiDose(ct,stf_singleModalityTwoJO,pln_2, cst, 0);
%%
stf_singleModalityTwo = matRad_generateStf(ct,cst,pln(2));
dij_singleModalityTwo = matRad_calcPhotonDose(ct,stf_singleModalityTwo,pln(2),cst)
%%
cst{6,6}{1} = DoseObjectives.matRad_SquaredDeviation(1,68);
cst{7,6}{1} = DoseObjectives.matRad_SquaredDeviation(2,56);
cst{1,6}{1} = DoseObjectives.matRad_EUDMin(3,30);
cst{8,6}{1} = DoseObjectives.matRad_EUDMin(4,30);
cst{2,6}{1} = DoseObjectives.matRad_MeanDose(5,30);
cst{4,6}{1} = DoseObjectives.matRad_EUDMin(6,30);
cst{9,6}{1} = DoseObjectives.matRad_MeanDose(7,20);
%%
test = matRad_
test = matRad_PriorityList1();


%%
cst{2,5}.Priority = 5;
cst{4,5}.Priority = 6;
cst{9,5}.Priority = 7;
%%

[resultGUI,resultGUIs,resultGUIs2,cstIt,cstIt2] = matRad_2pecOptimization(dij_singleModalityTwoJO,ct,cst,pln_2)
%%
for i = 1:numel(resultGUIs2)
    figure
    phDose = resultGUIs2{i}{1}.physicalDose;
    matRad_plotSliceWrapper(gca,ct,cst,1,phDose,3,round(pln_2.propStf.isoCenter(1,3)./ct.resolution.z),[],0.75,colorcube,[], [0 max([phDose(:)])],[]);
    title('')
end
%%

for i = 1:numel(resultGUIs)
    figure
    phDose = resultGUIs{i}{1}.physicalDose;
    matRad_plotSliceWrapper(gca,ct,cst,1,phDose,3,round(pln_2.propStf.isoCenter(1,3)./ct.resolution.z),[],0.75,colorcube,[], [0 max([phDose(:)])],[]);
    title('')
end
%%
for i = 1:numel(resultGUIs2)
    figure
    matRad_indicatorWrapperMixed(matRad_setOverlapPriorities(cst),pln_2,resultGUIs2{i})
end

