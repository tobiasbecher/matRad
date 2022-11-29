%% Example: Photon Treatment Plan
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this example we will show 
% (i) how to load patient data into matRad
% (ii) how to setup a photon dose calculation and 
% (iii) how to inversely optimize beamlet intensities
% (iv) how to visually and quantitatively evaluate the result

%% Patient Data Import
% Let's begin with a clear Matlab environment. Then, import the TG119
% phantom into your workspace. The phantom is comprised of a 'ct' and 'cst'
% structure defining the CT images and the structure set. Make sure the 
% matRad root directory with all its subdirectories is added to the Matlab 
% search path.

matRad_rc; %If this throws an error, run it from the parent directory first to set the paths
load('Liver.mat');

cst{15,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(0,0));
cst{16,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(100,45));
%%  
%cst2 = cst;
%cst{16,6}{1} = DoseObjectives.matRad_SquaredUnderdosing(800,45); 
%cst{16,6}{2} = DoseObjectives.matRad_SquaredOverdosing(800,48); 

%%
%size(cst2{16,6},2)
%% 
% The file TG119.mat contains two Matlab variables. Let's check what we 
% have just imported. First, the 'ct' variable comprises the ct cube along
%with some meta information describing properties of the ct cube (cube 
% dimensions, resolution, number of CT scenarios). Please note that 
%multiple ct cubes (e.g. 4D CT) can be stored in the cell array ct.cube{}
display(ct);
%%
% The 'cst' cell array defines volumes of interests along with information 
% required for optimization. Each row belongs to one certain volume of 
% interest (VOI), whereas each column defines different properties. 
% Specifically, the second and third column  show the name and the type of 
% the structure. The type can be set to OAR, TARGET or IGNORED. The fourth 
% column contains a linear index vector that lists all voxels belonging to 
% a certain VOI.
display(cst);
%%
% The fifth column represents meta parameters for optimization. The overlap
% priority is used to resolve ambiguities of overlapping structures (voxels 
% belonging to multiple structures will only be assigned to the VOI(s) with
% the highest overlap priority, i.e.. the lowest value). The parameters 
% alphaX and betaX correspond to the tissue's photon-radiosensitivity 
% parameter of the linear quadratic model. These parameter are required for
% biological treatment planning using a variable RBE. Let's output the meta 
% optimization parameter of the target, which is stored in the thrid row:

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% matlab structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this case
% we want to use photons. Then, we need to define a treatment machine to 
% correctly load the corresponding base data. matRad includes base data for
% generic photon linear accelerator called 'Generic'. By this means matRad 
% will look for 'photons_Generic.mat' in our root directory and will use 
% the data provided in there for dose calculation

pln.radiationMode = 'photons';  
pln.machine       = 'Generic';

%%
% Define the flavor of optimization along with the quantity that should be
% used for optimization. Possible quantities used for optimization are: 
% physicalDose: physical dose based optimization; 
% effect: biological effect based optimization;
% RBExD: RBE weighted dose based optimzation;
% Possible biological models are:
% none:        use no specific biological model
% constRBE:    use a constant RBE
% MCN:         use the variable RBE McNamara model for protons
% WED:         use the variable RBE Wedenberg model for protons
% LEM:         use the biophysical variable RBE Local Effect model for carbons
% As we are  using photons, we simply set the parameter to 'physicalDose' and
% and 'none'
quantityOpt    = 'physicalDose';                                     
modelName      = 'none';  

%%
% Now we have to set some beam parameters. We can define multiple beam 
% angles for the treatment and pass these to the plan as a vector. matRad 
% will then interpret the vector as multiple beams. In this case, we define
% linear spaced beams from 0 degree to 359 degree in 40 degree steps. This 
% results in 9 beams. All corresponding couch angles are set to 0 at this 
% point. Moreover, we set the bixelWidth to 5, which results in a beamlet 
% size of 5 x 5 mm in the isocenter plane. The number of fractions is set 
% to 30. Internally, matRad considers the fraction dose for optimization, 
% however, objetives and constraints are defined for the entire treatment.
pln.numOfFractions         = 30;
pln.propStf.gantryAngles   = [0:40:359];
pln.propStf.couchAngles    = zeros(1,numel(pln.propStf.gantryAngles));
pln.propStf.bixelWidth     = 5;

%%
% Obtain the number of beams and voxels from the existing variables and 
% calculate the iso-center which is per default the center of gravity of 
% all target voxels.
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%% dose calculation settings
% set resolution of dose calculation and optimization
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

%%
% Enable sequencing and disable direct aperture optimization (DAO) for now.
% A DAO optimization is shown in a seperate example.
pln.propOpt.runSequencing = 1;
pln.propOpt.runDAO        = 0;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen');

%%
% and et voila our treatment plan structure is ready. Lets have a look:
display(pln);


%% Generate Beam Geometry STF
% The steering file struct comprises the complete beam geometry along with 
% ray position, pencil beam positions and energies, source to axis distance (SAD) etc.
stf = matRad_generateStf(ct,cst,pln);

%%
% Let's display the beam geometry information of the 6th beam
display(stf(6));

%% Dose Calculation
% Let's generate dosimetric information by pre-computing dose influence 
% matrices for unit beamlet intensities. Having dose influences available 
% allows subsequent inverse optimization.
dij = matRad_calcPhotonDose(ct,stf,pln,cst);

%% Inverse Optimization for IMRT
% The goal of the fluence optimization is to find a set of beamlet/pencil 
% beam weights which yield the best possible dose distribution according to
% the clinical objectives and constraints underlying the radiation 
% treatment. Once the optimization has finished, trigger once the GUI to 
% visualize the optimized dose cubes.
%% Paretooptimization
% The goal of this step is to define a grid of penalty values that
% are then evaluated using the matRad_paretoGeneration method
% The VOI and their respective penalties are defined in the following way
% It is possible to have more than one objective function per VOI
% penVal stores the Grid which is then passed on. penGrid contains an
% version easier to visualize, however harder to loop over
%%
%[resultGUI2,o,wSt] = matRad_fluenceOptimization(dij,cst,pln);    
%%
VOI = {'Skin','PTV'};
nPoints = 20;
[pen,penGrid] = matRad_generateSphericalPenaltyGrid(nPoints,[1,1]);
%%
%matRad_plotPenaltyGrid(penGrid);

%%
penGrid = [[100,0];[0,100];penGrid]%;[0.9999,sqrt(1-0.9999^2)];[0.99999,sqrt(1-0.99999^2)]]*100;
%%
%%
%%
%objective function values are returned in order of ordering in VOI
nPoints = 22;
%returnStruct = matRad_paretoGeneration(dij,,pln,nPoints,VOI,[],penGrid);

returnStructCold = matRad_paretoGeneration(dij,cst,pln,nPoints,VOI,[],penGrid);
returnStructWarm = matRad_paretoGeneration(dij,cst,pln,nPoints,VOI,[],penGrid,true);
aaaaa
%%
returnStruct4.finds(4)
%%
[returnStruct3.weights(:,2),returnStruct2.weights(:,2)]
%%
resultGUI2.individualObj(2)
%%
returnStruct1.finds(2)
%%
returnStruct1.wInit(1)
%%
aaaaa
%%
returnStructWarmSq0.finds(6)
%%
matRad_plotParetoSurface(returnStructWarmSq0.finds,returnStructWarmSq0.penGrid/100,returnStructWarmSq0.VOIObj);

%%
aaaa
%%
load('resultsLiverSqOD0LinObjAll.mat')
load('resultsLiverSqOD0All.mat')
%%
penLin = returnStructAllLin0.penGrid;
penSq = returnStructAllSq0.penGrid;
findsLin = returnStructAllLin0.finds;
findsSq = returnStructAllSq0.finds;
VOIObjLin = returnStructAllLin0.VOIObj;
VOIObjSq = returnStructAllSq0.VOIObj;
%%
%load('resultsLiverPGENExtensive.mat')
%save('resultsLiverSqOD0LinObjAll.mat','-v7.3','returnStructAllLin0');
%load('resultsLiverSqOD0All.mat')
%load('resultsLiverSqOD0.mat');
%%
matRad_plotParetoSurface(returnStruct2.finds,returnStruct2.penGrid,returnStruct2.VOIObj);

%%
matRad_plotParetoSurface(returnStructAllSq0.finds,returnStructAllSq0.penGrid/100,returnStructAllSq0.VOIObj);

%%
matRad_plotParetoSurface(findsLin,penLin/100,VOIObjLin)
%%

matRad_plotParetoSurface(findsSq,penSq/100,VOIObjSq)
%%
figure
plot(findsLin(:,1).^2,findsLin(:,2).^2,'.',color = 'r')
hold on 
plot(findsSq(:,1),findsSq(:,2),'.',color = 'b')
%%
%load('resultsLiverPGEN2ObjectivesSqObj0D.mat')
%%
returnStruct2.finds
%%

matRad_plotParetoSurface([findsSq;returnStruct2.finds],[penSq/100;returnStruct2.penGrid],returnStruct2.VOIObj)





%%
resultGUI = matRad_calcCubes(returnStruct2.weights(:,2),dij);
resultGUI2 = matRad_calcCubes(returnStruct3.weights(:,2),dij);
%%
resultGUI3 = matRad_calcCubes(returnStruct4.weights(:,2),dij);
%%
%recalculate Plan from weights
resultGUIsOneSided = cell(size(weightsWarm,1),1);
for i =1:size(weightsWarm,2)
    w = weightsWarm(:,i);
    size(w)
    resultGUIsOneSided{i} = matRad_calcCubes(w,dij);
end

%%
resultGUIsWarm = cell(size(weightsCold,1),1);
for i =1:size(weightsCold,2)
    w = weightsWarm(:,i);
    size(w)
    resultGUIsWarm{i} = matRad_calcCubes(w,dij);
end
%%
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUIs{19});
%%
figure
imagesc(resultGUIsOneSided{1}.physicalDose(:,:,slice)),colorbar, colormap(jet);
%% Plot the Resulting Dose Slice

% Let's plot the transversal iso-center dose slice

slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
%%
for i=1:size(weightsWarm,2)
    figure
    imagesc(resultGUIsOneSided{i}.physicalDose(:,:,slice)),colorbar, colormap(jet);
end
%%
figure
'aaa'
plane      = 3;
absDiffCube = resultGUI3.physicalDose-resultGUI.physicalDose;
figure,title( 'fine beam spacing plan - coarse beam spacing plan')
matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],[],colorcube);