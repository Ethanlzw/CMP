clc;
clear;
close all;

%% Geometry and material
L = 1.0;
H = 2.0;
k = 0.6;
thickness = 1.0;

%% Mesh controls
nx = 16;
ny = 16;
nTriCols = floor(nx/2);   % 0 = pure Quad4, nx = pure Tri3

if nTriCols < 0 || nTriCols > nx
    error('nTriCols must satisfy 0 <= nTriCols <= nx.');
end

%% Fourier reference truncation
Nterms = 50;

%% Create FE model
model = struct();
model.fieldNames = {'T'};

model.elesType = {
    1, 'Tri3_Thermal',  1, 3;
    2, 'Quad4_Thermal', 1, 4
};

%% Generate nodes
dx = L / nx;
dy = H / ny;

numNodes = (nx + 1) * (ny + 1);
nodesInfo = zeros(numNodes, 3);

nodeID = 0;
for j = 0:ny
    for i = 0:nx
        nodeID = nodeID + 1;
        nodesInfo(nodeID,:) = [nodeID, i*dx, j*dy];
    end
end
model.nodesInfo = nodesInfo;

node = @(i,j) j*(nx+1) + i + 1;   % i=0..nx, j=0..ny

%% Generate mixed mesh
triElems = [];
quadElems = [];
eleID = 0;

for j = 0:(ny-1)
    for i = 0:(nx-1)
        n1 = node(i  , j  );
        n2 = node(i+1, j  );
        n3 = node(i+1, j+1);
        n4 = node(i  , j+1);

        if i < nTriCols
            eleID = eleID + 1;
            triElems = [triElems; eleID, n1, n2, n3];

            eleID = eleID + 1;
            triElems = [triElems; eleID, n1, n3, n4];
        else
            eleID = eleID + 1;
            quadElems = [quadElems; eleID, n1, n2, n3, n4];
        end
    end
end

model.elesInfo = struct();
model.elesInfo.elesSet1 = triElems;
model.elesInfo.elesSet2 = quadElems;
model.elesSetName = {'elesSet1','elesSet2'};

%% Problem data
model.problem = struct();
model.problem.thickness = thickness;

%% Section/material data
numEles = eleID;
model.section = zeros(numEles, 3);   % [eleID, kxx, kyy]
model.section(:,1) = 1:numEles;
for e = 1:numEles
    model.section(e,2:3) = [k, k];
end

%% Element loads
model.eleLoadData = cell(numEles,1);
for e = 1:numEles
    model.eleLoadData{e} = struct('heatSource', 0.0);
end

%% Global nodal vector
nNodeDoF = model.elesType{1,3};
model.nodalVector = zeros(size(model.nodesInfo,1)*nNodeDoF, 1);

%% Essential BCs

essBC = {}; % {[node list], [node temperature]}

% Left boundary x = 0
for j = 0:ny
    nid = node(0,j);
    essBC(end+1,1:2) = {nid, 0.0};
end

% Right boundary x = 1
for j = 0:ny
    nid = node(nx,j);
    essBC(end+1,1:2) = {nid, 0.1};
end

% Bottom boundary y = 0
% Avoid conflicting temperature at the left-bottom corner.
% Left boundary value is kept as the corner value.
for i = 1:nx
    nid = node(i,0);
    x = model.nodesInfo(nid,2);
    essBC(end+1,1:2) = {nid, 0.1};
end

% Top boundary y = H = 2
% Avoid conflicting temperature at the left-top corner.
% Left boundary value is kept as the corner value.
for i = 1:nx
    nid = node(i,ny);
    x = model.nodesInfo(nid,2);
    essBC(end+1,1:2) = {nid, 0.1};
end

model.essBC = essBC;

% Natural BCs
model.natBC = {}; % {[node list], [M,S]}

%% Solve FE problem
[globalMatrix, globalVector] = buildSystem(model);
bcData = buildBCData(model);
[reaction, solution] = solveSystem(globalMatrix, globalVector, bcData);

result.u = solution;
result.f = reaction;

%% Postprocessing
nodal = extractNodalFields(model, result.u);
results = postprocessModel(model, result.u);

%% Nodal averaged element fields
Tele = averageFieldToNodes(model, results, 'T');
qN = averageFieldToNodes(model, results, 'heatFlux');

%% Print results
disp('===== Fourier-Series Thermal Validation =====');
fprintf('nx = %d, ny = %d, nTriCols = %d\n', nx, ny, nTriCols);
fprintf('Number of nodes    = %d\n', size(model.nodesInfo,1));
fprintf('Number of elements = %d\n', numEles);
fprintf('Fourier terms      = %d\n', Nterms);

%% Plot
plotMesh(model, 'Fourier Thermal Validation Mesh', false, false);

plotNodalScalar(model, nodal.T, 'FE Temperature');

plotNodalScalar(model, Tele, 'Averaged Element Temperature');
plotNodalScalar(model, qN(:,1), 'Nodal Averaged q_x');
plotNodalScalar(model, qN(:,2), 'Nodal Averaged q_y');