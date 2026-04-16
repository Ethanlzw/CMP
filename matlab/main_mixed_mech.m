clc;
clear;
close all;

%% Geometry and material
L = 10.0;
H = 1.0;
t = 1.0;

E  = 210e9;
nu = 0.30;
P  = 1000.0;

%% Mesh controls
nx = 20;
ny = 4;
nTriCols = nx;

if nTriCols < 0 || nTriCols > nx
    error('nTriCols must satisfy 0 <= nTriCols <= nx.');
end

%% Create FE model
model = struct();
model.fieldNames = {'ux','uy'};

model.elesType = {
    1, 'Tri3_Mech',  2, 3;
    2, 'Quad4_Mech', 2, 4
};

model.elesSetName = {'setTri', 'setQuad'};

%% Material / problem / loads
model.material = [E, nu];
model.problem = struct();
model.problem.thickness = t;
model.problem.mechType = 'planeStress';

model.loads = struct();
model.loads.bodyForce = [0; 0];

%% Generate nodes and elements
% TODO: fill nodesInfo and elesInfo
model.nodesInfo = [];
model.elesInfo = struct();
model.elesInfo.setTri = [];
model.elesInfo.setQuad = [];

%% Boundary conditions definition
model.bcDef = struct();

% Essential BC example
model.bcDef.essential = struct();
model.bcDef.essential.nodes  = [];
model.bcDef.essential.fields = {};
model.bcDef.essential.values = [];

% Natural BC example
model.bcDef.natural = struct();
model.bcDef.natural.edges  = [];
model.bcDef.natural.type   = {};
model.bcDef.natural.values = [];
model.bcDef.natural.setName = {};

%% Build BC data
bcData = buildBCData(model);
model.bcData = bcData;

%% Build and solve
[globalMatrix, globalVector] = buildSystem(model);
[u, reaction] = solveSystem(globalMatrix, bcData, globalVector);

%% Post-process
result = postprocessModel(model, u);

%% Plot
% plotMesh(model);
% plotDeformedMesh(model, result, 1.0);
% plotDeformedMeshContour(model, result, 'uy', 1.0);