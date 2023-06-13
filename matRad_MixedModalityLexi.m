matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 50000;
  matRad_cfg.propOpt.defaultAccChangeTol = 1e-7;
%%
%{
load 'TG119.mat'

%%
cst{1,6}{1} = struct(DoseObjectives.matRad_MeanDose(100,25));
cst{3,6}{1} = struct(DoseObjectives.matRad_MeanDose(100,30));
cst{2,6}{1} = struct(DoseConstraints.matRad_MinMaxDose(45,57.5));
%%
%}
load('TG119_MOD.mat')
%}
%%
%%
% 
% meta information for treatment plan (1) 
pln(1).numOfFractions  = 5;
pln(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln(1).machine         = 'Generic';

% beam geometry settings
pln(1).propStf.bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
pln(1).propStf.gantryAngles    = [0]; % [?] ;
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
pln(1).propDoseCalc.doseGrid.resolution.x = 8; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.y = 8; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.z = 8; % [mm]
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
% 
% meta information for treatment plan (2) 

pln(2).numOfFractions  = 25;
pln(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln(2).machine         = 'Generic';

% beam geometry settings
pln(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(2).propStf.gantryAngles    = [0:70:359]; % [?] ;
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
pln(2).propDoseCalc.doseGrid.resolution.x = 8; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.y = 8; % [mm]   
pln(2).propDoseCalc.doseGrid.resolution.z = 8; % [mm]
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

cst = matRad_prepCst(cst, sparecst);
% Plan Wrapper
plnJO = matRad_plnWrapper(pln);
% Stf Wrapper
stf = matRad_stfWrapper(ct,cst,plnJO);
% Dij Calculation
%%
dij = matRad_calcCombiDose(ct,stf,plnJO,cst,false);
% Fluence optimization 
%%


%%
load('TG119.mat')
cst{1,6}{1} = DoseObjectives.matRad_MeanDose(100,0);
%cst{2,6}{1} = DoseConstraints.matRad_MinMaxDose(45,57.5);
%cst{2,6}{1}.epsilon = 3e-2;
cst{3,6} = [];
%%
resultGUI = matRad_fluenceOptimizationJO(dij,cst,plnJO);

physicalDose = resultGUI{1}.physicalDose*5 + resultGUI{2}.physicalDose*25;
%%


%%
figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(resultGUI{1}.physicalDose(:,:,slice)*5)%-physicalDose2(:,:,slice))
colorbar()
    
figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(resultGUI{2}.physicalDose(:,:,slice)*25)%-physicalDose2(:,:,slice))
colorbar()

figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(physicalDose(:,:,slice))
colorbar()

%%
wOpt = [resultGUI{1}.wUnsequenced;resultGUI{2}.wUnsequenced];
%%
delta = 1.03;
load(['TG119_MOD.mat'])
%%
const1 = DoseConstraints.matRad_ObjectiveConstraint(cst{1,6}{1},25*delta,0);

cst{1,6}{1} = const1;
%%
resultGUI2 = matRad_fluenceOptimizationJO(dij,cst,plnJO,wOpt);
%% plots

physicalDose = resultGUI2{1}.physicalDose*5 + resultGUI2{2}.physicalDose*25;
figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(resultGUI2{1}.physicalDose(:,:,slice)*5)%-physicalDose2(:,:,slice))
colorbar()
    
figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(resultGUI2{2}.physicalDose(:,:,slice)*25)%-physicalDose2(:,:,slice))
colorbar()

figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(physicalDose(:,:,slice))
colorbar()
%%
wOpt = [resultGUI2{1}.wUnsequenced;resultGUI2{2}.wUnsequenced]

load('TG119.mat')
const1 = DoseConstraints.matRad_ObjectiveConstraint(cst{1,6}{1},resultGUI2{3}.objectives*1.03,0);
%%
cst{1,6}{1} = const1;
cst{2,6}{1} = const2;
%%
resultGUI3 = matRad_fluenceOptimizationJO(dij,cst,plnJO,wOpt);
%% plots

physicalDose = resultGUI3{1}.physicalDose*5 + resultGUI3{2}.physicalDose*25;

figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(resultGUI3{1}.physicalDose(:,:,slice)*5)%-physicalDose2(:,:,slice))
colorbar()
    
figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(resultGUI3{2}.physicalDose(:,:,slice)*25)%-physicalDose2(:,:,slice))
colorbar()

figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(physicalDose(:,:,slice))
colorbar()
%%
