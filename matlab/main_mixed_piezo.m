clc;
clear;
close all;

%% Geometry and material
L = 1.0;
H = 1.0;
t = 1.0;

E   = 210e9;
nu  = 0.30;
eps = 1.0e-8;
e31 = -5.0;

%% Model
model = struct();
model.fieldNames = {'ux','uy','phi'};

model.elesType = {
    1, 'Tri3_Piezo',  3, 3;
    2, 'Quad4_Piezo', 3, 4
};

model.elesSetName = {'setTri', 'setQuad'};

model.material = [E, nu, eps, e31];
model.problem = struct();
model.problem.thickness = t;
model.problem.mechType = 'planeStress';

model.loads = struct();
model.loads.bodyForce = [0; 0];
model.loads.chargeSource = 0.0;

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