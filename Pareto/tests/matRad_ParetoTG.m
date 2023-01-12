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
%global timed;
%timed = [];
%% Patient Data Import
% Let's begin with a clear Matlab environment. Then, import the TG119
% phantom into your workspace. The phantom is comprised of a 'ct' and 'cst'
% structure defining the CT images and the structure set. Make sure the 
% matRad root directory with all its subdirectories is added to the Matlab 
% search path.

matRad_rc; %If this throws an error, run it from the parent directory first to set the paths
%%
load('TG119.mat');
%%
%cst{2,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(45,55));
cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(0,0));
cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(0,0));
%cst{1,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(0,38));
%cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(0,43));
%cst{2,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(45,55));
%cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(0,0));

pln.radiationMode = 'photons';  
pln.machine       = 'Generic';

quantityOpt    = 'physicalDose';                                     
modelName      = 'none';  


pln.numOfFractions         = 30;
pln.propStf.gantryAngles   = [0:40:359];
pln.propStf.couchAngles    = zeros(1,numel(pln.propStf.gantryAngles));
pln.propStf.bixelWidth     = 5;


pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

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
VOI = {'Core','OuterTarget','BODY'};
%%
%objective function values are returned in order of ordering in VOI
%returnStruct = matRad_paretoGenerationPGEN(dij,cst,pln,VOI);
%

nPoints = 100;
[pen,penGrid] = matRad_generateSphericalPenaltyGrid(nPoints,[1,1,1]);
%%

penGrid = [[100/sqrt(3),100/sqrt(3),100/sqrt(3)];[100,0,0];[0,0,100];[0,100,0];penGrid];
matRad_plotPenaltyGrid(penGrid);
nPoints = size(penGrid,1);
%%
%returnStruct = matRad_paretoGenerationPGEN(dij,cst,pln,VOI)
%cst{1,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(0,38));
%cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(0,43));
cst{1,6} = [];
cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(0,0));
cst{1,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(0,38));

cst{2,6} = [];
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(1000,50));
%cst{2,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(45,55));

cst{3,6} = [];
cst{3,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(0,0));
cst{3,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(0,43))

returnStructCold = matRad_paretoGeneration(dij,cst,pln,nPoints,VOI,[],penGrid);


%returnStruct = matRad_paretoGenerationPGEN(dij,cst,pln,VOI);
aaaaaaaa
%%
i = 4;
[k,facets,normals,cs,dists,w]= matRad_convexHull(returnStruct.finds,returnStruct.penGrid)
%%
P = returnStruct.finds;
[k,vol] = convhulln(P);
trisurf(k,P(:,1),P(:,2),P(:,3),'FaceColor','cyan')
%%
%save('resultsTG1191constraint.mat','-v7.3','returnStructCold');
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
i = 4;
matRad_plotParetoSurface(returnStructCold.finds(1:i,:),returnStructCold.penGrid(1:i,:)/100,VOI)
%%


%returnStruct2.penGrid
returnStructCold.penGrid(1:4,:)

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