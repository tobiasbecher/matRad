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
load('Prostate.mat');

%%
cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,0));
cst{8,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,0));
cst{6,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(60,78));
cst{1,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(0,60));
cst{8,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(0,60));
cst{9,6}{1} = struct(DoseConstraints.matRad_MinMaxDose(0,50));
cst{7,6}{1} = struct(DoseConstraints.matRad_MinMaxDose(52,66));
%%
%% 
pln.radiationMode = 'photons';  
pln.machine       = 'Generic';

%%
quantityOpt    = 'physicalDose';                                     
modelName      = 'none';  
pln.numOfFractions         = 30;
pln.propStf.gantryAngles   = [0:40:359];
pln.propStf.couchAngles    = zeros(1,numel(pln.propStf.gantryAngles));
pln.propStf.bixelWidth     = 5;

pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

pln.propDoseCalc.doseGrid.resolution.x = 6; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 6; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 6; % [mm]

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


%%
matRadGUI;
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
VOI = {'Rectum','PTV_68','Bladder'};
%%

%cst{6,6} = [];
%cst{6,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(1000,50));
%cst{6,6}{2} = struct(DoseConstraints.matRad_MinMaxDose(62,78));

%penGrid = [[100/sqrt(3),100/sqrt(3),100/sqrt(3)];[100,0,0]]%;[0,100,0];[0,0,100]];%;penGrid];

%nPoints = size(penGrid,1);


%%
returnStruct = matRad_RennenRadio(dij,cst,pln,VOI);
%returnStruct2 = matRad_paretoGenerationPGEN(dij,cst,pln,VOI);
aaaaaaaa

%%
i = 5;
matRad_plotParetoSurface(returnStruct.finds(1:i,:),returnStruct.penGrid(1:i,:),VOI)
%%
returnStruct.finds
%%
ps = returnStruct.finds

[k,facets] = matRad_ParetoSurfFromFacets(ps)
figure
trisurf(k,ps(:,1),ps(:,2),ps(:,3),'FaceColor','cyan')

figure
trisurf(facets(all(facets,2),:),ps(:,1),ps(:,2),ps(:,3),'FaceColor','cyan')
hold on 
scatter3(ps(:,1),ps(:,2),ps(:,3),'MarkerEdgeColor','black',...
        'MarkerFaceColor',[0 0 0])


%%

weights = returnStruct.weights;
resultGUIs = cell(size(weights,4));
for i =1:5
    w = weights(:,i);
    resultGUIs{i} = matRad_calcCubes(w,dij);
end
%%
resultGUIs{4}
%%

slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
for i=1:5
    figure
    imagesc(resultGUIs{i}.physicalDose(:,:,slice)),colorbar, colormap(jet);
end
%%
resultGUI = resultGUIs{2};
matRadGUI;
%%
penGrid(1,:)
penGrid(50,:)
%%
'aaa'
plane      = 3;
absDiffCube = resultGUIs{31}.physicalDose-resultGUIs{1}.physicalDose;
figure,title( 'fine beam spacing plan - coarse beam spacing plan')
matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],[],colorcube);