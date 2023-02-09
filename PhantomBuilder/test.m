clear all;
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
'a' 
objective = struct(DoseObjectives.matRad_SquaredDeviation(800,60));
testobj = CubicObj('Volume1','OAR',objective,[40,40,40],[0,0,0]);  
testobj2 = SphericalObj('Volume2','TARGET',objective,30,[0,0,0]);
testobj3 = SphericalObj('Volume3','OAR',objective,40,[0,0,0]);  
%%
cst = {}
%%
'a'
cst = testobj.initializeParameters(ct,cst)  
cst = testobj2.initializeParameters(ct,cst)
cst = testobj3.initializeParameters(ct,cst)
%%

%%
for i = 1:size(cst,1)
    vIxVOI = cst{i,4}{1};
    ct.cubeHU{1}(vIxVOI) = 0; % assign HU of water
end
%%

matRadGUI


