clc;
clear;
close all;

%% Geometry and material
plateLength = 1.0;
plateHeight = 1.0;
plateThickness = 1.0;
thermalConductivity = 10.0;

%% Mesh controls
nX = 20;
nY = 10;
nTriColumns = nX;

if nTriColumns < 0 || nTriColumns > nX
    error('nTriColumns must satisfy 0 <= nTriColumns <= nX.');
end

%% Create model
model = struct();
model.fieldNames = {'T'};

model.elesType = {
    1, 'Tri3_Thermal',  1, 3;
    2, 'Quad4_Thermal', 1, 4
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

%% Section
nElements = elementID;
model.section = zeros(nElements, 2);
model.section(:,1) = 1:nElements;
model.section(:,2) = thermalConductivity;

%% Element loads
model.eleLoadData = cell(nElements,1);
for e = 1:nElements
    model.eleLoadData{e} = struct('heatSource', 0.0);
end

model.nodalVector = zeros(size(model.nodesInfo,1), 1);

%% Essential BC
leftNodeIDs = [];
rightNodeIDs = [];
for j = 0:nY
    leftNodeIDs = [leftNodeIDs, node(0,j)];
    rightNodeIDs = [rightNodeIDs, node(nX,j)];
end

model.essBC = {
    leftNodeIDs,  [100];
    rightNodeIDs, [0]
};

%% No natural BC
model.natBC = {};

%% Build and solve
[globalMatrix, globalVector] = buildSystem(model);
bcData = buildBCData(model);
[reaction, solution] = solveSystem(globalMatrix, globalVector, bcData);

result.u = solution;
result.f = reaction;

nodal = extractNodalFields(model, result.u);
results = postprocessModel(model, result.u);
gradTAtNodes = averageFieldToNodes(model, results, 'tempGrad');

plotMesh(model, 'Mixed Thermal Mesh', false, false);
plotNodalScalar(model, nodal.T, 'Temperature');
plotNodalScalar(model, sqrt(gradTAtNodes(:,1).^2 + gradTAtNodes(:,2).^2), 'Temperature Gradient Magnitude');