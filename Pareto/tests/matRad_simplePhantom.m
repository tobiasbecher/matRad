
%% set matRad runtime configuration
clear all; %somewhat needed for the phantom builder
matRad_rc; %If this throws an error, run it from the parent directory first to set the paths

ctDim = [200,200,100]; % y,x,z dimensions
ctResolution = [2,2,3]; % y,x,z the same here!

builder = matRad_PhantomBuilder(ctDim,ctResolution,1)

%define objectives for the VOI

objective1 = struct(DoseObjectives.matRad_SquaredDeviation(800,45));
objective2 = struct(DoseObjectives.matRad_SquaredOverdosing(400,0));
objective3 = struct(DoseObjectives.matRad_SquaredOverdosing(10,0));     

builder.addSphericalTarget('Volume1',20,'objectives',objective1,'HU',0);
builder.addBoxOAR('Volume2',[60,30,60],'offset',[0 -15 0],'objectives',objective2);
builder.addBoxOAR('Volume3',[60,30,60],'offset',[0 15 0],'objectives',objective3);

%% Get the ct and cst (stored as properties of the phantomBuilder class)

[ct,cst] = builder.getctcst();
%%
%%

%matRadGUI;
%%

%%
pln.radiationMode = 'protons';            
pln.machine       = 'Generic';

%%
% Define the biological optimization model for treatment planning along
% with the quantity that should be used for optimization. Possible model values 
% are:
% 'none':     physical optimization;
% 'constRBE': constant RBE of 1.1; 
% 'MCN':      McNamara-variable RBE model for protons; 
% 'WED':      Wedenberg-variable RBE model for protons
% 'LEM':      Local Effect Model 
% and possible quantityOpt are 'physicalDose', 'effect' or 'RBExD'.
modelName    = 'none';
quantityOpt  = 'physicalDose';                                             

%%
% The remaining plan parameters are set like in the previous example files
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [0 180];
pln.propStf.couchAngles   = zeros(size(pln.propStf.gantryAngles));
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
%pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1)*[200 200 75]
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve nominal scenario for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen'); 

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);
%%


cst{1,6}{2} = DoseConstraints.matRad_MinMaxDose(42.75,51.7)

%%
resultGUI = matRad_fluenceOptimization(dij,cst,pln)
matRadGUI
%%
returnStructTG100 = matRad_RennenRadio(dij,cst,pln,100);
%%
load('dataFull.mat')
data.finds
returnStructCold.finds
%%

matRad_sliderUI(data,dij,pln,ct)
%%
aaaa
%%
data = returnStructTG100
%%
%Visualize surface
%data = retStruct
ps =  data.finds
[k2,vol] = convhulln(ps)
[k,facets] = matRad_ParetoSurfFromFacets(ps)
figure
trisurf(k2,ps(:,1),ps(:,2),ps(:,3),'FaceColor','red')
hold on
trisurf(facets(all(facets,2),:),ps(:,1),ps(:,2),ps(:,3),'FaceColor','blue')
scatter3(data.finds(:,1),data.finds(:,2),data.finds(:,3))
%%
%matRad_plotPenaltyGrid(returnStructCold.penGrid)


%% Recalculate plan from stored weights

resultGUIs ={}
for i = 1:size(returnStructCold.weights,2)
    resultGUIs{i} = matRad_calcCubes(returnStructCold.weights(:,i),dij);
end
%%
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
for i= 1:10
    figure
    imagesc(resultGUIs{i}.physicalDose(:,:,slice)),colorbar, colormap(jet);
    title(returnStructCold.penGrid(i,:));
end
%%
resultGUI = resultGUIs{3}
matRadGUI
%%