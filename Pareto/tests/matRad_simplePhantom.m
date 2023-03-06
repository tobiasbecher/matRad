clear all;
matRad_rc;
xDim = 200;
yDim = 200;
zDim = 100;

ct.cubeDim      = [yDim xDim zDim]; % second cube dimension represents the x-coordinate
ct.resolution.x = 2;
ct.resolution.y = 2;
ct.resolution.z = 3;
ct.numOfCtScen  = 1;
 
% create a ct image series with zeros - it will be filled later
ct.cubeHU{1} = ones(ct.cubeDim) * -1000;
%%


objective3 = struct(DoseObjectives.matRad_SquaredDeviation(800,45));
objective2 = struct(DoseObjectives.matRad_SquaredOverdosing(10,0));
objective1 = struct(DoseObjectives.matRad_SquaredOverdosing(400,0));

%%
builder = matRad_PhantomBuilder(ct);
builder.addCubicOAR('Volume1',[60,30,60],'offset',[0 -15 0],'objectives',objective1);
builder.addCubicOAR('Volume2',[60,30,60],'offset',[0 15 0],'objectives',objective2);
builder.addSphericalTarget('Volume3',20,'objectives',objective3)
%%
builder.updatecst()
%%  
builder.updatect()
%%
cst = builder.cst;
ct = builder.ct;
%
%%
cst{1,5}
%%
cst{1,5}.Priority = 3;
cst{2,5}.Priority = 2;
cst{3,5}.Priority = 1,

cst{3,6}{2} = DoseConstraints.matRad_MinMaxDose(42.75,51.7);
%{
cst{1,6}{2} = DoseConstraints.matRad_MinMaxDose(0,40);

cst{2,6}{2} = DoseConstraints.matRad_MinMaxDose(0,40);
%}
%%
%%

%matRadGUI;
%%

%%
pln.radiationMode = 'photons';            
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
pln.propStf.gantryAngles  = [0:70:355];
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
pln.propDoseCalc.doseGrid.resolution.x = 4; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 4; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 4; % [mm]

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
dij = matRad_calcPhotonDose(ct,stf,pln,cst);
%%
%{
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot the resulting dose slice
plane      = 3;
slice      = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
doseWindow = [0 max([resultGUI.physicalDose(:)])];

figure,title('phantom plan')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],[],colorcube,[],doseWindow,[]);


matRadGUI
%}

%%



%%
VOI = {'Volume1','Volume2','Volume3'};

%returnStructCold = matRad_RennenRadio(dij,cst,pln,VOI);
%%
%size(returnStructCold.weights(:,1))

%resultGUIs = {};
%size(returnStructCold.weights)


%%

%%
%data = returnStructCold;
%%

%%
load("data.mat")
matRad_plotParetoSurface(data.finds,data.penGrid,{'a','b','c'})


%%
load('dataFull.mat')

aaaa
%%

ps =  data.finds
[k2,vol] = convhulln(ps)
[k,facets,allnormals] = matRad_ParetoSurfFromFacets(ps)
figure
trisurf(k2,ps(:,1),ps(:,2),ps(:,3),'FaceColor','red')
hold on
trisurf(facets(all(facets,2),:),ps(:,1),ps(:,2),ps(:,3),'FaceColor','blue')
scatter3(data.finds(:,1),data.finds(:,2),data.finds(:,3))
%%
%matRad_plotPenaltyGrid(returnStructCold.penGrid)


L = min(ps,[],1);
U = max(ps,[],1);
ps = (ps-L)./(U-L);
ps
%%
matRad_plotPenaltyGrid(returnStructCold.penGrid)
%%
matRad_plotParetoSurface(returnStructCold.finds,returnStructCold.penGrid,VOI)
%%
ps = returnStructCold.finds;
[k,vol] = convhulln(ps);
trisurf(k,ps(:,1),ps(:,2),ps(:,3),'FaceColor','cyan')
%%
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
matRad_plotParetoSurface(data.finds,data.penGrid,{'a','b','c'})
%%
matRad_plotPenaltyGrid(data.penGrid)

%%
fS = data.finds(1:2,:)
%%
ls = data.removed{1}(end,:)
%%
ls
%%
aa = [fS;ls];

[ha,haha,no] = matRad_normalFromFacet(aa,1:3,1)
no
%%
zw = no.*(U'-L')
zw = zw/sqrt(sum(zw.^2))
%%
L = min(aa,[],1);
U = max(aa,[],1);
%%
aa = (aa-L)./(U-L);
%%
aa
[fS;ls]
%%
[ha,haha,no] = matRad_normalFromFacet(aa,1:3,1)