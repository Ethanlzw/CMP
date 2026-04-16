clc;
clear;
close all;

%% Geometry and material
L = 1.0;
H = 1.0;
k = 10.0;

%% Mesh controls
nx = 10;
ny = 10;
nTriCols = nx;

%% Create FE model
model = struct();
model.fieldNames = {'T'};

model.elesType = {
    1, 'Tri3_Thermal',  1, 3;
    2, 'Quad4_Thermal', 1, 4
};

model.elesSetName = {'setTri', 'setQuad'};

%% Material / problem / loads
model.material = k;
model.problem = struct();
model.problem.thickness = 1.0;

model.loads = struct();
model.loads.source = 0.0;

%% Mesh
model.nodesInfo = [];
model.elesInfo = struct();
model.elesInfo.setTri = [];
model.elesInfo.setQuad = [];

%% BC
model.bcDef = struct();
model.bcDef.essential = struct();
model.bcDef.essential.nodes  = [];
model.bcDef.essential.fields = {};
model.bcDef.essential.values = [];

model.bcDef.natural = struct();
model.bcDef.natural.edges  = [];
model.bcDef.natural.type   = {};
model.bcDef.natural.values = [];
model.bcDef.natural.setName = {};

%% Solve
bcData = buildBCData(model);
model.bcData = bcData;

[globalMatrix, globalVector] = buildSystem(model);
[u, reaction] = solveSystem(globalMatrix, bcData, globalVector);

result = postprocessModel(model, u);

%% Plot
% plotNodalScalar(model, result, 'T');