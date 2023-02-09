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
%%
load('Liver.mat');
%%
%cst{15,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,0));
%cst{16,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(1000,45));

%cst{15,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(0,45));
cst{16,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(40,50));
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
VOI = {'Skin','PTV'};
%%
%objective function values are returned in order of ordering in VOI
%returnStruct = matRad_paretoGenerationPGEN(dij,cst,pln,VOI);
%
returnStruct2 = matRad_paretoGenerationPGEN(dij,cst,pln,VOI);
aaaaaaaa
%%
save('resultsLiverPGEN1constr.mat','-v7.3','returnStruct2');
%%
returnStruct2
returnStructPGEN = returnStruct2;
%%
returnStruct2.finds
%%
returnStruct2.penGrid
%%
convhulln(returnStruct2.finds(1:i,:))
%%
matRad_plotParetoSurface(returnStruct2.finds(1:i,:),returnStruct2.penGrid(1:i,:),VOI)

%%

i = 5
[a,b,c,d,e,f] = matRad_convexHull(returnStruct2.finds(1:i,:),returnStruct2.penGrid(1:i,:));
%%
%returnStruct2.penGrid
returnStruct2.finds

%%
returnStructAllSq0.penGrid
%%
aaaaaaaaaaaaaa
%%
returnStructLiverPGEN2LinObj0 = returnStruct2
%save('resultsLiverPGEN2ObjectivesLinObj0D.mat','-v7.3','returnStructLiverPGEN2LinObj0');
%%
load('resultsLiverSqOD0All.mat')
%%

matRad_plotParetoSurface(returnStructAllSq0.finds,returnStructAllSq0.penGrid/100,returnStructAllSq0.VOIObj)
%%
matRad_plotParetoSurface(returnStruct2.finds,returnStruct2.penGrid,returnStruct2.VOIObj)
%%

returnStruct2.finds
returnStruct2.penGrid
%%
%%
penGrid = returnStruct2.penGrid;
fInd = returnStruct2.finds;
VOIObj = returnStruct2.VOIObj;
%%
matRad_plotParetoSurface(fInd,penGrid,VOIObj)
%%
n = 3;
fVals = fInd(1:n,:);
penVals = penGrid(1:n,:);
[k,facets,normals,cs,dists,w2] = matRad_convexHull(fVals,penVals);
penVals
w2
%%

n = 7;
fVals = fInd(1:n,:);
penVals = penGrid(1:n,:);
[k,facets,normals,cs,dists,w2] = matRad_convexHull(fVals,penVals);
penVals
w2
%%

w2
%%
penVals
%%
%for i = 1:size(k,1)
    plot(fVals(k',1),fVals(k',2),'color','k')
%end
%%
returnStruct2.finds
returnStruct2.findsDirect
%%
fVals = returnStruct2.finds;
weights = returnStruct2.penGrid;
%%
size(fVals)
size(weights)
%%
[k,normals,cs,dists,w] = matRad_convexHull(fVals,weights)
%%
k
%%
matRad_plotParetoSurface(fVals,penGrid,returnStruct2.VOIObj);
hold on
test = fVals(k,:)
plot(test(:,1),test(:,2),'-*')
%%


slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
for i=1:2:50
    figure
    imagesc(resultGUIs{i}.physicalDose(:,:,slice)),colorbar, colormap(jet);
end
%%
%%
penGrid(1,:)
penGrid(50,:)
%%
'aaa'
plane      = 3;
absDiffCube = resultGUIs{31}.physicalDose-resultGUIs{1}.physicalDose;
figure,title( 'fine beam spacing plan - coarse beam spacing plan')
matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],[],colorcube);