%% Load patient data
matRad_rc;
matRad_cfg = MatRad_Config.instance();
load 'TG119_MOD.mat'
%% Define meta info for pln(1)

%%
n = 30
pln(1).numOfFractions  = n;
pln(1).radiationMode   = 'carbon';           % either photons / protons / helium / carbon
pln(1).machine         = 'Generic';

% beam geometry settings
pln(1).propStf.bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
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
pln(2).numOfFractions  = n;
pln(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln(2).machine         = 'Generic';

% beam geometry settings
pln(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(2).propStf.gantryAngles    = [0:70:359];%[0:90:359]; % [?] ;
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

%cst = matRad_updatecst(cst,[0,0,1])

%%

%% Single Modality 1 JO
%%
%{
cst{1,6}{1}.parameters{1} = cst{1,6}{1}.parameters{1}/(30/n);
cst{2,6}{1}.parameters{1} = cst{2,6}{1}.parameters{1}/(30/n);
cst{3,6}{1}.parameters{1} = cst{3,6}{1}.parameters{1}/(30/n);
%}
%{
%cst{2,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(42,57.5));
cst{2,6}{2} = struct(DoseConstraints.matRad_MinMaxDVH(40,0.85,1))
%pln_1.scaling = 4;
%%



stf_singleModalityOne =  matRad_generateStf(ct,cst,pln(1),0);
dij_singleModalityOne = matRad_calcParticleDose(ct,stf_singleModalityOne,pln(1), cst, 0);
%%
resultGUI_singleModalityOne = matRad_fluenceOptimization(dij_singleModalityOne, cst,pln(1));

    %%
pln_1 = matRad_plnWrapper(pln(1));

% Stf
stf_singleModalityOneJO = matRad_stfWrapper(ct,cst, pln_1);

% dij
dij_singleModalityOneJO = matRad_calcCombiDose(ct,stf_singleModalityOneJO, pln_1, cst, 0);

% fluence optimization
%
resultGUI_singleModalityOneJO = matRad_fluenceOptimizationJO(dij_singleModalityOneJO, cst, pln_1);
%%
%}
%% Single Modality 2 JO

pln_2 = matRad_plnWrapper(pln(2));

% Stf
stf_singleModalityTwoJO = matRad_stfWrapper(ct,cst, pln_2);

% dij
dij_singleModalityTwoJO = matRad_calcCombiDose(ct,stf_singleModalityTwoJO,pln_2, cst, 0);
%%
cst{1,6}{1} = DoseObjectives.matRad_MeanDose(100,0);
cst{2,6}{1} = DoseConstraints.matRad_MinMaxDose(45,57.5);
cst{2,6}{1}.epsilon = 3e-2;
cst{3,6} = [];
%%
% fluence optimization
resultGUI_singleModalityTwoJO = matRad_fluenceOptimizationJO(dij_singleModalityTwoJO, cst, pln_2,wOpt);
'a'

%%
load('TG119_MOD.mat')
cst{1,6}{1} = DoseConstraints.matRad_MinMaxMeanDose(0,resultGUI_singleModalityTwoJO{1,2}.objectives*1.03);
cst{2,6}{1} = DoseConstraints.matRad_MinMaxDose(45,57.5);
cst{2,6}{1}.epsilon = 3e-2;
cst{3,6}{1} = DoseObjectives.matRad_MeanDose(100,0);
%%
wOpt = resultGUI_singleModalityTwoJO{1}.w
%%
resultGUI_singleModalityTwoJOTwo = matRad_fluenceOptimizationJO(dij_singleModalityTwoJO, cst, pln_2,wOpt);
%% Joint optimization
% 
% plnJO = matRad_plnWrapper(pln);
% 
% % Stf Wrapper
% stf = matRad_stfWrapper(ct,cst,plnJO);
% 
% % Dij Calculation
% dij = matRad_calcCombiDose(ct,stf,plnJO,cst,false);
% 
% % Fluence optimization 
% resultGUI = matRad_fluenceOptimizationJO(dij,cst,plnJO)
%%
mean(physicalDose2(cst{1,4}{1}))*30

%% Native single modality 1


%% Native signle modality 2
stf_singleModalityTwo =  matRad_generateStf(ct,cst,pln(2),0);
dij_singleModalityTwo = matRad_calcPhotonDose(ct,stf_singleModalityTwo,pln(2), cst, 0);
resultGUI_singleModalityTwo = matRad_fluenceOptimization(dij_singleModalityTwo, cst,pln(2));

%%

physicalDose =  resultGUI_singleModalityOneJO{1}.physicalDose; %+ resultGUI_singleModalityOneJO{2}.physicalDose*25;
%%
physicalDose2 =  resultGUI_singleModalityTwoJO{1}.physicalDose;%*5 + resultGUI_singleModalityTwoJO{2}.physicalDose*25;
%%

physicalDose22 =  resultGUI_singleModalityTwoJOTwo{1}.physicalDose;%*5 + resultGUI_singleModalityTwoJO{2}.physicalDose*25;
%%
physicalDose3 =  resultGUI_singleModalityOne.physicalDose;%*5 + resultGUI_singleModalityTwoJO{2}.physicalDose*25;
%%
physicalDose4 =  resultGUI_singleModalityTwo.physicalDose;%*5 + resultGUI_singleModalityTwoJO{2}.physicalDose*25;
%physicalDose2 = resultGUINew{1}.physicalDose*5 + resultGUINew{2}.physicalDose*25;
%%
figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(physicalDose(:,:,slice))%-physicalDose2(:,:,slice))
colorbar()
%%
figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(30*physicalDose2(:,:,slice))%-physicalDose2(:,:,slice))
colorbar()
%%
figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(30*physicalDose22(:,:,slice))%-physicalDose2(:,:,slice))
colorbar()
%%
figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(physicalDose3(:,:,slice))%-physicalDose2(:,:,slice))
colorbar()
%%
figure
slice = round(pln(1).propStf.isoCenter(1,3)./ct.resolution.z);
imagesc(physicalDose4(:,:,slice))%-physicalDose2(:,:,slice))
colorbar()



%% Comparison
slice = round(size(resultGUI_singleModalityOne.physicalDose,3)/2);
profile = 70;
cMap = matRad_getColormap('gammaIndex');
f = figure;
subplot(2,2,1);
    imagesc(resultGUI_singleModalityOne.physicalDose(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title('Single modality')
subplot(2,2,2);
    imagesc(resultGUI_singleModalityOneJO{1}.physicalDose(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title('JO');
subplot(2,2,3);
    absDifference = zeros(size(resultGUI_singleModalityOne.physicalDose, [1,2]));
    absDifference = 100*(resultGUI_singleModalityOne.physicalDose(:,:,slice) - resultGUI_singleModalityOneJO{1}.physicalDose(:,:,slice))./resultGUI_singleModalityOne.physicalDose(:,:,slice);
    absDifference(resultGUI_singleModalityOne.physicalDose(:,:,slice)== 0) = 0;
    imagesc(absDifference);
    clim([0 100]);
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title('% difference');
    colorbar;

ax = subplot(2,2,4);
    gamma_SingleModalityOne = matRad_gammaIndex(resultGUI_singleModalityOne.physicalDose, resultGUI_singleModalityOneJO{1}.physicalDose, [3 3 2.5], [3 3]);
    passingVoxels = squeeze(sum(gamma_SingleModalityOne(:,:,slice) < 1, [1,2]));
    passingRate_SingleModalityOne = 100*passingVoxels./numel(gamma_SingleModalityOne(:,:,slice));

    imagesc(gamma_SingleModalityOne(:,:,slice));
    clim([0 2]);
    colormap(ax,cMap);
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);

    title(['Slice ', num2str(slice), ', passing rate = ', num2str(passingRate_SingleModalityOne)]);
    colorbar;
sgtitle(['Proton physical dose, slice = ', num2str(slice)]);



f = figure;
subplot(2,2,1);
    imagesc(resultGUI_singleModalityTwo.physicalDose(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title('Single modality')
subplot(2,2,2);
    imagesc(resultGUI_singleModalityTwoJO{1}.physicalDose(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title('JO');
subplot(2,2,3);
    absDifference = zeros(size(resultGUI_singleModalityTwo.physicalDose, [1,2]));
    absDifference = 100*(resultGUI_singleModalityTwo.physicalDose(:,:,slice) - resultGUI_singleModalityTwoJO{1}.physicalDose(:,:,slice))./resultGUI_singleModalityTwo.physicalDose(:,:,slice);
    absDifference(resultGUI_singleModalityTwo.physicalDose(:,:,slice)== 0) = 0;
    imagesc(absDifference);
    clim([0 100]);
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title('% difference');
    colorbar;

ax = subplot(2,2,4);
    gamma_SingleModalityTwo = matRad_gammaIndex(resultGUI_singleModalityTwo.physicalDose, resultGUI_singleModalityTwoJO{1}.physicalDose, [3 3 2.5], [3 3]);
    passingVoxels = squeeze(sum(gamma_SingleModalityTwo(:,:,slice) < 1, [1,2]));
    passingRate_SingleModalityTwo = 100*passingVoxels./numel(gamma_SingleModalityTwo(:,:,slice));

    imagesc(gamma_SingleModalityTwo(:,:,slice));
    clim([0 2]);
    colormap(ax,cMap);
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);

    title(['Slice ', num2str(slice), ', passing rate = ', num2str(passingRate_SingleModalityTwo)]);
    colorbar;
sgtitle(['Photon physical dose, slice = ', num2str(slice)]);
% f = figure;
% subplot(1,2,1);
% imagesc(resultGUI_singleModalityTwo.physicalDose(:,:,slice));
% matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
% title('Single modality')
% subplot(1,2,2);
% 
% imagesc(resultGUI_singleModalityTwoJO{1}.physicalDose(:,:,slice));
% matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
% title('JO');
% sgtitle(['Photon physical dose, slice = ', num2str(slice)]);
% %colorbar;

%% gamma analysis
slice_gamma = [50, 60, 70, 80];
gamma_SingleModalityOne = matRad_gammaIndex(resultGUI_singleModalityOne.physicalDose, resultGUI_singleModalityOneJO{1}.physicalDose, [3 3 2.5], [3 3]);
passingVoxels = squeeze(sum(gamma_SingleModalityOne(:,:,slice_gamma) < 1, [1,2]));
passingRate_SingleModalityOne = 100*passingVoxels./numel(gamma_SingleModalityOne(:,:,slice_gamma(1)));

f = figure;
for k=1:size(slice_gamma,2)
    subplot(2,2,k)
    imagesc(gamma_SingleModalityOne(:,:,slice_gamma(k)));
    clim([0 2]);
    colormap(cMap);
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice_gamma(k));
    title(['Slice ', num2str(slice_gamma(k)), ', passing rate = ', num2str(passingRate_SingleModalityOne(k))]);
    colorbar;
end

sgtitle(['Gamma index Protons']);


gamma_SingleModalityTwo = matRad_gammaIndex(resultGUI_singleModalityTwo.physicalDose, resultGUI_singleModalityTwoJO{1}.physicalDose, [3 3 2.5], [3 3]);
passingVoxels = squeeze(sum(gamma_SingleModalityTwo(:,:,slice_gamma) < 1, [1,2]));
passingRate_SingleModalityTwo = 100*passingVoxels/numel(gamma_SingleModalityTwo(:,:,slice_gamma(1)));

f = figure;
for k=1:size(slice_gamma,2)
subplot(2,2,k)
    imagesc(gamma_SingleModalityTwo(:,:,slice_gamma(k)));
    clim([0 2]);
    colormap(cMap);
    title(['Slice ', num2str(slice_gamma(k)), ', passing rate = ', num2str(passingRate_SingleModalityTwo(k))]);
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice_gamma(k));
    colorbar;
end
sgtitle(['Gamma index Photons']);
%}