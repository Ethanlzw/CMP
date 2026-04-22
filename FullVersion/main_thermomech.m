clc;
clear;
close all;

%% Geometry and material
plateLength = 1.0;
plateHeight = 0.2;
plateThickness = 1.0;

E = 210e9;
nu = 0.30;
thermalConductivity = 10.0;
thermalExpansion = 1.2e-5;

%% Mesh controls
nX = 20;
nY = 4;
nTriColumns = nX;

if nTriColumns < 0 || nTriColumns > nX
    error('nTriColumns must satisfy 0 <= nTriColumns <= nX.');
end

%% Create model
model = struct();
model.fieldNames = {'ux','uy','T'};

model.elesType = {
    1, 'Tri3_ThermoMech',  3, 3;
    2, 'Quad4_ThermoMech', 3, 4
};

%% Generate nodes
dx = plateLength / nX;
dy = plateHeight / nY;

nNodes = (nX + 1) * (nY + 1);
nodesInfo = zeros(nNodes, 3);

nodeID = 0;
for j = 0:nY
    for i = 0:nX
        nodeID = nodeID + 1;
        nodesInfo(nodeID,:) = [nodeID, i*dx, j*dy];
    end
end
model.nodesInfo = nodesInfo;

node = @(i,j) j*(nX+1) + i + 1;

%% Mesh
triElements = [];
quadElements = [];
elementID = 0;

for j = 0:(nY-1)
    for i = 0:(nX-1)
        n1 = node(i  , j  );
        n2 = node(i+1, j  );
        n3 = node(i+1, j+1);
        n4 = node(i  , j+1);

        if i < nTriColumns
            elementID = elementID + 1;
            triElements = [triElements; elementID, n1, n2, n3];

            elementID = elementID + 1;
            triElements = [triElements; elementID, n1, n3, n4];
        else
            elementID = elementID + 1;
            quadElements = [quadElements; elementID, n1, n2, n3, n4];
        end
    end
end

model.elesInfo = struct();
model.elesInfo.elesSet1 = triElements;
model.elesInfo.elesSet2 = quadElements;
model.elesSetName = fieldnames(model.elesInfo);

%% Problem data
model.problem = struct();
model.problem.thickness = plateThickness;
model.problem.mechType = 'planeStress';
model.problem.deltaTRef = 0.0;

%% Section
nElements = elementID;
model.section = zeros(nElements, 5);
model.section(:,1) = 1:nElements;

for e = 1:nElements
    model.section(e,2:5) = [E, nu, thermalConductivity, thermalExpansion];
end

%% Element loads
model.eleLoadData = cell(nElements,1);
for e = 1:nElements
    model.eleLoadData{e} = struct('bodyForce', [0; 0], 'heatSource', 0.0);
end

model.nodalVector = zeros(size(model.nodesInfo,1)*3, 1);

%% Essential BC
leftNodeIDs = [];
rightNodeIDs = [];
for j = 0:nY
    leftNodeIDs = [leftNodeIDs, node(0,j)];
    rightNodeIDs = [rightNodeIDs, node(nX,j)];
end

% Left edge fixed in mechanics and prescribed temperature 100
% Right edge temperature 0
model.essBC = {
    leftNodeIDs,  [0, 0, 100];
    rightNodeIDs, [NaN, NaN, 0]
};

%% No natural BC
model.natBC = {};

%% Build and solve
[globalMatrix, globalVector] = buildSystem(model);
bcData = buildBCData(model);
[reaction, solution] = solveSystem(globalMatrix, globalVector, bcData);

result.u = solution;
result.f = reaction;

%% Postprocessing
nodal = extractNodalFields(model, result.u);
results = postprocessModel(model, result.u);
vmAtNodes = averageFieldToNodes(model, results, 'vonMises');

%% Plot
plotMesh(model, 'Thermo-Mechanical Mesh', false, false);
plotNodalScalar(model, nodal.T, 'Temperature');
plotNodalScalar(model, vmAtNodes, 'Von Mises Stress');
plotDeformedMesh(model, result.u, 1e4, 'Thermo-Mechanical Deformed Mesh');
plotDeformedMeshContour(model, result.u, 1e4, 'Thermo-Mechanical Displacement Contour');