clear all;
matRad_rc;
xDim = 30;
yDim = 50;
zDim = 30;
%%
matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 50000;
matRad_cfg.propOpt.defaultAccChangeTol = 1e-7;
%%
ct.cubeDim      = [yDim xDim zDim]; % second cube dimension represents the x-coordinate
ct.resolution.x = 2;
ct.resolution.y = 3;
ct.resolution.z = 3;
ct.numOfCtScen  = 1;
 
% create a ct image series with zeros - it will be filled later
ct.cubeHU{1} = ones(ct.cubeDim) * -1000;%%
%ct.cubeHU{1} = zeros(ct.cubeDim)
objective3 = struct(DoseObjectives.matRad_SquaredDeviation(800,45));
objective2 = struct(DoseObjectives.matRad_SquaredOverdosing(10,0));
objective1 = struct(DoseObjectives.matRad_SquaredOverdosing(400,0));
%%
builder = matRad_PhantomBuilder(ct.cubeDim,[2,3,3],1);
builder.addBoxOAR('Volume1',[14,20,20],'offset',[-14,0,0],'objectives',objective1);
builder.addBoxOAR('Volume2',[10,20,20],'offset',[10,0,0],'objectives',objective2);
builder.addBoxTarget('Volume3',[10,10,20],'offset',[-1,0,0],'objectives',objective3);
builder.addBoxOAR('Volume4',[40,26,20],'offset',[-3,0,0])
%builder.addSphericalTarget('Volume3',10,'offset',5,'objectives',objective3)
%%

%%
[ct,cst] = builder.getctcst();
%%
%%
cst{1,5}.Priority = 3;
cst{2,5}.Priority = 2;
cst{3,5}.Priority = 1;
%%
%%%%
pln(1).numOfFractions  = 5;
pln(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln(1).machine         = 'Generic';

% beam geometry settings
pln(1).propStf.bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
pln(1).propStf.gantryAngles    = [0];
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
pln(2).propStf.gantryAngles    = [0:80:359];%[0:90:359]; % [?] ;
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
%%
%[resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln);
%%

pln_2 = matRad_plnWrapper(pln(2));

% Stf
stf_singleModalityTwoJO = matRad_stfWrapper(ct,cst, pln_2);

% dij
dij_singleModalityTwoJO = matRad_calcCombiDose(ct,stf_singleModalityTwoJO,pln_2, cst, 0);
%%
%resultGUI2 = matRad_fluenceOptimizationJO(dij_singleModalityTwoJO,cst,pln_2)
%%
 pln_M= matRad_plnWrapper(pln);

% Stf
stf_M = matRad_stfWrapper(ct,cst, pln_M);

% dij
dij_M = matRad_calcCombiDose(ct,stf_M, pln_M, cst, 0);
%%
cst{3,6}{1} = DoseConstraints.matRad_MinMaxDose(45,57.5);
cst{1,6}{1} = DoseObjectives.matRad_MeanDose(2,10);
cst{2,6}{1} = DoseObjectives.matRad_MeanDose(1,10);
%%
stf_singleModalityTwo =  matRad_generateStf(ct,cst,pln(2),0);
dij_singleModalityTwo = matRad_calcPhotonDose(ct,stf_singleModalityTwo,pln(2), cst, 0);
%%
[resultGUI,resultGUIs,resultGUIs2,cstL] = matRad_2pecOptimization2(dij_singleModalityTwo,ct,cst,pln(2));

%%
[resultGUI,resultGUIs,resultGUIs2,cstL] = matRad_2pecOptimizationMixed(dij_M,cst,pln_M);
%{
resultGUI = matRad_fluenceOptimizationJO(dij_M,cst,pln_M)
%%

physicalDose2 = resultGUI{1}.physicalDose*5 + resultGUI{2}.physicalDose*25;
%%

figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(resultGUI{1}.physicalDose(:,:,slice)*5)%-physicalDose2(:,:,slice))
colorbar()
    
%%
figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(resultGUI{2}.physicalDose(:,:,slice)*25)%-physicalDose2(:,:,slice))
colorbar()
    %%
figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(physicalDose2(:,:,slice))%-physicalDose2(:,:,slice))
colorbar()
%}
%% conversion for GUI

%resultGUI = resultGUI2{1}
pln2 = pln
pln = pln2(2)
matRadGUI
