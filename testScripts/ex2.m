
matRad_rc; %If this throws an error, run it from the parent directory first to set the paths
load('TG119.mat');

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

ixTarget = 3;
display(cst{ixTarget,5});

%%
% The sixth column contains optimization information such as objectives and
% constraints which are required to calculate the objective function value. 
% Please note, that multiple objectives/constraints can be defined for
% individual structures. Here, we have defined a squared deviation 
% objective making it 'expensive/costly' for the optimizer to over- and 
% underdose the target structure (both are equally important). 

display(cst{ixTarget,6});

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
pln.propStf.gantryAngles   = [0:70:359];
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
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

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
% display(stf(6));
%% Dose Calculation
% Let's generate dosimetric information by pre-computing dose influence 
% matrices for unit beamlet intensities. Having dose influences available 
% allows subsequent inverse optimization.
dij = matRad_calcPhotonDose(ct,stf,pln,cst);
%%
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
%%
matRadGUI