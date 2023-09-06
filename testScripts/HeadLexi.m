%% Load patient data
matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 50000;
matRad_cfg.propOpt.defaultAccChangeTol = 1e-7;
%%
load 'HEAD_AND_NECK.mat'
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

%% Prepare cst

cst = matRad_prepCst(cst, sparecst);

%cst = matRad_updatecst(cst,[0,0,1])

%%
%{
pln_M= matRad_plnWrapper(pln);

% Stf
stf_M = matRad_stfWrapper(ct,cst, pln_M);

% dij
dij_M = matRad_calcCombiDose(ct,stf_M, pln_M, cst, 0);
%}
%%
%{
cst{15,6}{1} = DoseConstraints.matRad_MinMaxDose(55,77);
cst{16,6}{1} = DoseConstraints.matRad_MinMaxDose(60,82);
%cst{15,6}{2} = DoseObjectives.matRad_SquaredOverdosing(100,72.45);

%cst{16,6}{2} = DoseObjectives.matRad_SquaredUnderdosing(100,63);
cst{13,6}{1} = DoseObjectives.matRad_EUDMin(1,26,4);
%cst{13,6}{1} = DoseObjectives.matRad_MeanDose(1,20);
cst{14,6}{1} = DoseObjectives.matRad_EUDMin(2,26,4);
%cst{14,6} = [];
cst{17,6}{1} = DoseObjectives.matRad_MeanDose(3,28);
%}
%[resultGUI,resultGUIs,resultGUIs2,cstL] = matRad_2pecOptimization(dij_M,cst,pln_M,wInit);
%%
%resultGUI = matRad_fluenceOptimizationJO(dij_M,cst,pln_M)
%%
%wInit = [resultGUI{1}.w;resultGUI{2}.w];
%%


%resultGUI = matRad_fluenceOptimizationJO(dij_M,cst,pln_M,wInit)
%%
%resultGUI = matRad_fluenceOptimizationJO(dij_singleModalityTwoJO,cstIt2,pln_2)




%% Single Modality 1 JO
%%

%%
%{


stf_singleModalityOne =  matRad_generateStf(ct,cst,pln(1),0);
dij_singleModalityOne = matRad_calcParticleDose(ct,stf_singleModalityOne,pln(1), cst, 0);
%%

cst{15,6}{1} = DoseConstraints.matRad_MinMaxDose(50,85);
cst{16,6}{1} = DoseConstraints.matRad_MinMaxDose(50,85);
%cst{15,6}{2} = DoseObjectives.matRad_SquaredOverdosing(100,72.45);
%cst{16,6}{2} = DoseObjectives.matRad_SquaredUnderdosing(100,63);
cst{13,6}{1} = DoseObjectives.matRad_MeanDose(100,0);
%cst{14,6}{1} = DoseObjectives.matRad_MeanDose(100,40);
cst{14,6} = [];
cst{17,6} = [];

resultGUI_singleModalityOne = matRad_fluenceOptimization(dij_singleModalityOne, cst,pln(1));

%%
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
figure
imagesc(resultGUI_singleModalityOne.physicalDose(:,:,slice))
colorbar()

%%
%}
%{
pln_1 = matRad_plnWrapper(pln(1));

% Stf
stf_singleModalityOneJO = matRad_stfWrapper(ct,cst, pln_1);

% dij
dij_singleModalityOneJO = matRad_calcCombiDose(ct,stf_singleModalityOneJO, pln_1, cst, 0);

cst{15,6}{1} = DoseConstraints.matRad_MinMaxDose(55,77);
cst{16,6}{1} = DoseConstraints.matRad_MinMaxDose(60,82);
%cst{15,6}{2} = DoseObjectives.matRad_SquaredOverdosing(100,72.45);

%cst{16,6}{2} = DoseObjectives.matRad_SquaredUnderdosing(100,63);
cst{13,6}{1} = DoseObjectives.matRad_EUDMin(1,26,4);
%cst{13,6}{1} = DoseObjectives.matRad_MeanDose(1,20);
cst{14,6}{1} = DoseObjectives.matRad_EUDMin(2,26,4);
%cst{14,6} = [];
cst{17,6}{1} = DoseObjectives.matRad_MeanDose(3,28);
%%
%resultGUI = matRad_fluenceOptimization(dij,cst,pln)

[resultGUI,resultGUIs,resultGUIs2,cstL] = matRad_2pecOptimization(dij_singleModalityOneJO,cst,pln_1);

%%


for i = 1:numel(resultGUIs2)
    i
    physicalDose =  resultGUIs2{i}{1}.physicalDose;%*5 + resultGUI_singleModalityTwoJO{2}.physicalDose*25;
    figure
    title('i')
    slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
    imagesc(physicalDose(:,:,slice))%-physicalDose2(:,:,slice))
    colorbar()
    figure
    dvh = matRad_calcDVH(cst,physicalDose)
    matRad_showDVH(dvh,cst,pln_1);
    
end
%}
%%

pln_2 = matRad_plnWrapper(pln(2));

% Stf
stf_singleModalityTwoJO = matRad_stfWrapper(ct,cst, pln_2);

% dij
dij_singleModalityTwoJO = matRad_calcCombiDose(ct,stf_singleModalityTwoJO,pln_2, cst, 0);
%%


'a'
%%

%stf_singleModalityTwo =  matRad_generateStf(ct,cst,pln(2),0);
%dij_singleModalityTwo = matRad_calcPhotonDose(ct,stf_singleModalityTwo,pln(2), cst, 0);
%%
cst{15,6}{1} = DoseConstraints.matRad_MinMaxDose(55,77);
cst{16,6}{1} = DoseConstraints.matRad_MinMaxDose(60,82);
%cst{15,6}{2} = DoseObjectives.matRad_SquaredOverdosing(100,72.45);

%cst{16,6}{2} = DoseObjectives.matRad_SquaredUnderdosing(100,63);
cst{13,6}{1} = DoseObjectives.matRad_EUDMin(1,26,4);
%cst{13,6}{1} = DoseObjectives.matRad_MeanDose(1,20);
cst{14,6}{1} = DoseObjectives.matRad_EUDMin(2,26,4);
%cst{14,6} = [];
cst{17,6}{1} = DoseObjectives.matRad_MeanDose(3,28);
%%
%[resultGUI,resultGUIs,resultGUIs2,cstL] = matRad_2pecOptimization(dij_singleModalityTwoJO,cst,pln_2);
load('HeadNeckFirstTryEUDPhotonsMixed.mat')
%%
matRad_indicatorWrapperMixed(matRad_setOverlapPriorities(cstL),pln_2,resultGUIs2{2})

%%
resultGUIMod = resultGUIs2{2};
%resultGUIMod = resultGUI{1}
resultGUIMod.physicalDose = resultGUIMod.physicalDose*30
matRad_indicatorWrapper(matRad_setOverlapPriorities(cst),pln_2,resultGUIMod)
%%
matRad_indicatorWrapperMixed(matRad_setOverlapPriorities(cst),pln_2,resultGUIs2{1})
%%
matRad_indicatorWrapperMixed(matRad_setOverlapPriorities(cst),pln_2,resultGUIs2{2})

%%
testGUI = resultGUIs2{2}
%%
d = testGUI{1}.physicalDose;
%%
cstO = matRad_setOverlapPriorities(cst)
%%
dM = d(cstO{14,4}{1})*30;
%%
n = 3.5
1/(numel(dM))^(1/n) * sum(dM.^n)^(1/n)
%%
mean(dM)
%%
for i = 1:numel(resultGUIs2)
    physicalDose =  resultGUIs2{i}{1}.physicalDose;%*5 + resultGUI_singleModalityTwoJO{2}.physicalDose*25;
    figure
    title('i')
    slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
    imagesc(physicalDose(:,:,slice))%-physicalDose2(:,:,slice))
    colorbar()
    figure
    dvh = matRad_calcDVH(cst,physicalDose);
    matRad_showDVH(dvh,cst,pln_2);
end