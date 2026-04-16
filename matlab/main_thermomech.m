clc;
clear;
close all;

%% Geometry and material
L = 1.0;
H = 1.0;
t = 1.0;

E     = 210e9;
nu    = 0.30;
k     = 20.0;
alpha = 1.2e-5;

%% Model
model = struct();
model.fieldNames = {'ux','uy','T'};

model.elesType = {
    1, 'Tri3_ThermoMech',  3, 3;
    2, 'Quad4_ThermoMech', 3, 4
};

model.elesSetName = {'setTri', 'setQuad'};

model.material = [E, nu, k, alpha];
model.problem = struct();
model.problem.thickness = t;
model.problem.mechType = 'planeStress';
model.problem.Tref = 0.0;

model.loads = struct();
model.loads.bodyForce = [0; 0];
model.loads.source = 0.0;

model.nodesInfo = [];
model.elesInfo = struct();
model.elesInfo.setTri = [];
model.elesInfo.setQuad = [];

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

bcData = buildBCData(model);
model.bcData = bcData;

[globalMatrix, globalVector] = buildSystem(model);
[u, reaction] = solveSystem(globalMatrix, bcData, globalVector);

result = postprocessModel(model, u);