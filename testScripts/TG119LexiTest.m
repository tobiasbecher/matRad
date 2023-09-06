
matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 50000;
matRad_cfg.propOpt.defaultAccChangeTol = 1e-2;
%%
load('TG119.mat');
%%

pln(1).numOfFractions  = 5;
pln(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln(1).machine         = 'Generic';

% beam geometry settings
pln(1).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(1).propStf.gantryAngles    = [45 315];
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
pln(2).propStf.bixelWidth      = 5 %; % [mm] / also corresponds to lateral spot spacing for particles
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
pln(2).propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.y = 3; % [mm]
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
%% Generate Beam Geometry STF

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
dij_M = matRad_calcCombiDose(ct,stf_M, pln_M, cst, 0);
%%

%%

test = matRad_PriorityList1();
test.addObjective(1,DoseObjectives.matRad_SquaredDeviation(100,50),4,2);
test.addObjective(2,DoseObjectives.matRad_EUDMin(100,0,5),25,1);
test.addObjective(3,DoseObjectives.matRad_MeanDose(100,0),15,3);
test.addObjective(4,DoseObjectives.matRad_SquaredOverdosing(100,30),0,3);



cst{2,6}{1} = DoseObjectives.matRad_SquaredDeviation(1,50);
cst{1,6}{1} = DoseObjectives.matRad_EUDMin(2,25,5);

%cst{3,6}{2} = DoseObjectives.matRad_SquaredOverdosing(3,50);
cst{3,6}{1} = DoseObjectives.matRad_MeanDose(3,15);
cst{3,6}{2} = DoseObjectives.matRad_SquaredOverdosing(4,30);


%cst{2,6}{1} = DoseConstraints.matRad_MinMaxDose(45,57.5);
%cst{1,6}{2} = DoseConstraints.matRad_MinMaxDose(0,45);
%cst{3,6}{2} = DoseConstraints.matRad_MinMaxDose(0,45);
'a'
[resultGUIs,resultGUIs2,cst1,cst2]=  matRad_2pecOptimizationNew(test,dij_singleModalityTwoJO,ct,cst,pln_2);   

[resultGUIP,resultGUIsP,resultGUIs2P,cstLP] = matRad_2pecOptimization(dij_singleModalityTwoJO,ct,cst,pln_2);
%%
%[resultGUI,resultGUIs,resultGUIs2,cstL] = matRad_2pecOptimizationMixed(dij_M,cst,pln_M);
%%
%load('TG1193ObjNice.mat')
%%
%load('TG119MixedModsLexi.mat')
%%
for i = 1:numel(resultGUIs2)
    matRad_indicatorWrapperMixed(cst,pln_M,resultGUIs2{i})
end
%%
for i = 1:numel(resultGUIs2P)
    phDose2 = resultGUIs2P{i}{1}.physicalDose;
    max(phDose(:))
end
%%
matRad_plotLexicographicMixed(pln_M,ct,cst,resultGUIs2)
%%
%%
for i = 1:numel(resultGUIs2)
    figure
    phDose = resultGUIs2{i}{1}.physicalDose*5+resultGUIs2{i}{2}.physicalDose*25;
    matRad_plotSliceWrapper(gca,ct,cst,1,phDose,3,round(pln_2.propStf.isoCenter(1,3)./ct.resolution.z),[],0.75,colorcube,[], [0 max([phDose(:)])],[]);
    title('')
end
%%
phDose = resultGUIs2P{1}{1}.physicalDose;
phDose3 = resultGUIs2P{3}{1}.physicalDose;
matRad_plotSliceWrapper(gca,ct,cst,1,phDose3-phDose,3,round(pln_2.propStf.isoCenter(1,3)./ct.resolution.z),[],0.75,colorcube,[], [min([phDose3(:)-phDose(:)]) max([phDose3(:)-phDose(:)])],[]);
  
%%