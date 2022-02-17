%% Example: Proton LET optimization
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% In this example we will show 


%% Patient Data Import
% Let's begin with a clear Matlab environment and import the prostate
% patient into your workspace

matRad_rc; %If this throws an error, run it from the parent directory first to set the paths
clear all;
load('TG119.mat');

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines the most 
% important cornerstones of your treatment plan.

%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this
% example we would like to use protons for treatment planning. Next, we
% need to define a treatment machine to correctly load the corresponding 
% base data. matRad features generic base data in the file
% 'proton_Generic.mat'; consequently the machine has to be set accordingly
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

%%
% We use no biological optimization (to combine physical dose, XBD(LET) and
% XBD(DADR) sensibly
pln.propOpt.bioOptimization = 'none';

  
%%
% Now we have to set the remaining plan parameters.
pln.numOfFractions        = 1;
pln.propStf.gantryAngles  = [180];
pln.propStf.couchAngles   = [0];
pln.propStf.bixelWidth    = 5;
pln.propStf.longitudinalSpotSpacing = 3;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
% for particles with tabulated LET it is possible to also calculate the LET 
% disutribution alongside the physical dose. Therefore you need to activate
% the corresponding option during dose calculcation
pln.propDoseCalc.calcLET = true;

% Lets generate dosimetric information by pre-computing dose influence 
% matrices for unit beamlet intensities. Having dose influences available 
% allows for subsequent inverse optimization. 
dij = matRad_calcParticleDose(ct,stf,pln,cst);

%Set a fixed current in the dij
dij.fixedCurrent = 5;

%Set some MU data
dij.minMU = 0; % minimum MU to be produced
dij.MU = 5; %factor MU -> w (now means 1 MU = 0.16 * w = 0.16*1e6 particles)

%If you want to output RBE-weighted dose as well:
dij.RBE = 1.1;

%% Here you can modify values used in the computation of XBD(DADR) and
% XBD(LET);
%dij.c = 0.04; %FOr XBD(LET)
%dij.xbd_DADR_k = 0.5;
%dij.xbd_DADR_t = 40;
%dij.xbd_DADR_a_num = 4; %a = a_num / DADR_t

%% Initial optimization without DADR objectives
cst{1,6}{1}.parameters{1} = 4;
cst{2,6}{1}.parameters{1} = 10;
cst{3,6}{1}.parameters{1} = 5;

resultGUI = matRad_fluenceOptimization(dij,cst,pln);
resultGUI = matRad_calcDADR(dij,resultGUI);

%% Plot dose and DADR
plotF = figure;
pDose = subplot(2,3,1);
matRad_plotSliceWrapper(pDose,ct,cst,1,resultGUI.physicalDose,3,65);

pDADR = subplot(2,3,2);
matRad_plotSliceWrapper(pDADR,ct,cst,1,resultGUI.DADR,3,65); %unit should be Gy/s

pLET = subplot(2,3,3);
matRad_plotSliceWrapper(pLET,ct,cst,1,resultGUI.LET,3,65);

dvh_opt1 = matRad_calcDVH(cst,resultGUI.DADR);


%% Run DADR optimization
%Add Dose-Rate objective on Core
XBD_DADR_ref = 0.5*cst{1,6}{1}.parameters{1}; %Half of the dose prescribed above
cst{1,6}{2} = struct(XBDDADRObjectives.matRad_SquaredUnderXBDDADR(10000, XBD_DADR_ref));

XBD_LET_ref = cst{2,6}{1}.parameters{1} * 1.1; %GO towards an RBE of 1.1
cst{2,6}{2} = struct(XBDLETObjectives.matRad_SquaredDeviationXBDLET(5000, XBD_LET_ref));

%RUn optimization
resultGUI = matRad_fluenceOptimization(dij,cst,pln,resultGUI.w);
resultGUI = matRad_calcDADR(dij,resultGUI);

%% Plot dose and DADR
figure(plotF);
pDose = subplot(2,3,4);
matRad_plotSliceWrapper(pDose,ct,cst,1,resultGUI.physicalDose,3,65);

pDADR = subplot(2,3,5);
matRad_plotSliceWrapper(pDADR,ct,cst,1,resultGUI.DADR,3,65); %unit should be Gy/s

pLET = subplot(2,3,6);
matRad_plotSliceWrapper(pLET,ct,cst,1,resultGUI.LET,3,65);

%% DADRs
dvh_opt2 = matRad_calcDVH(cst,resultGUI.DADR);
figure;
matRad_showDVH(dvh_opt1,cst,pln); hold on;
matRad_showDVH(dvh_opt2,cst,pln,2);

%% Inverse Optimization for IMPT
matRadGUI




