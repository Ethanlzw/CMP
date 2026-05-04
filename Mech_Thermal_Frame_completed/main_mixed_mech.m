clc;
clear;
close all;

%% Geometry and material
L = 10.0;          % beam length
H = 1.0;           % beam height
t = 1.0;           % thickness

E  = 210e9;
nu = 0.30;
P  = 1000.0;       % total downward load at right edge

%% Mesh controls
nx = 80;            % number of cells along x
ny = 16;            % number of cells along y
nTriCols = nx;      % left part uses triangles, right part uses quads

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

% Helper for node numbering
node = @(i,j) j*(nx+1) + i + 1;   % i=0..nx, j=0..ny

%% Generate element connectivity
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
            % split one quad cell into two triangles
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

model.elesSetName = fieldnames(model.elesInfo);

%% Problem data
model.problem = struct();
model.problem.thickness = t;
model.problem.mechType = 'planeStress';

%% Section/material data
numEles = eleID;
model.section = zeros(numEles, 3);   % [eleID, E, nu]
model.section(:,1) = 1:numEles;
for e = 1:numEles
    model.section(e,2:3) = [E, nu];
end

%% Element loads
model.eleLoadData = cell(numEles,1);
for e = 1:numEles
    model.eleLoadData{e} = struct('bodyForce', [0; 0]);
end

%% Global nodal vector
nNodeDoF = model.elesType{1,3};
model.nodalVector = zeros(size(model.nodesInfo,1)*nNodeDoF, 1);

%% Boundary conditions
% Fix left edge
leftNodes = [];
for j = 0:ny
    leftNodes = [leftNodes, node(0,j)];
end

model.essBC = {
    leftNodes, [0,0]
};

% Uniform downward traction on right edge so total resultant = P
rightNodes = [];
for j = 0:ny
    rightNodes = [rightNodes, node(nx,j)];
end

ty = -P / (H * t);
model.natBC = {
    rightNodes, [0, ty]
};

%% Build and solve
[globalMatrix, globalVector] = buildSystem(model);
bcData = buildBCData(model);
[reaction, solution] = solveSystem(globalMatrix, globalVector, bcData);

result.u = solution;
result.f = reaction;

%% Postprocessing
nodal = extractNodalFields(model, result.u);
results = postprocessModel(model, result.u);

stressN = averageFieldToNodes(model, results, 'stress');
vmN = averageFieldToNodes(model, results, 'vonMises');

%% Beam theory references
A = t * H;
I = t * H^3 / 12;
G = E / (2 * (1 + nu));
kappa = 5/6;

w_tip_EB  = P * L^3 / (3 * E * I);
w_tip_Tim = P * L^3 / (3 * E * I) + P * L / (kappa * G * A);

sigma_top_exact = -P * L * (H/2) / I;
sigma_bot_exact = +P * L * (H/2) / I;

%% FE quantities of interest
% Use the mid-height node on the right edge if available, otherwise average all right-edge uy
rightNodeIDs = rightNodes(:);
rightCoordsY = model.nodesInfo(rightNodeIDs, 3);

[~, midIdx] = min(abs(rightCoordsY - H/2));
midRightNode = rightNodeIDs(midIdx);

uy_tip_mid = nodal.uy(midRightNode);
uy_tip_avg = mean(nodal.uy(rightNodeIDs));

% Total vertical reaction on left edge
fixed_dofs_y = zeros(numel(leftNodes),1);
for k = 1:numel(leftNodes)
    fixed_dofs_y(k) = (leftNodes(k)-1)*2 + 2;
end
Ry_total = sum(result.f(fixed_dofs_y));

%% Errors
err_tip_EB_mid  = abs(abs(uy_tip_mid) - w_tip_EB)  / abs(w_tip_EB);
err_tip_Tim_mid = abs(abs(uy_tip_mid) - w_tip_Tim) / abs(w_tip_Tim);

err_tip_EB_avg  = abs(abs(uy_tip_avg) - w_tip_EB)  / abs(w_tip_EB);
err_tip_Tim_avg = abs(abs(uy_tip_avg) - w_tip_Tim) / abs(w_tip_Tim);

err_reaction = abs(abs(Ry_total) - P) / P;

%% Print results
disp('===== Refined Cantilever Beam Validation =====');
fprintf('nx = %d, ny = %d, nTriCols = %d\n', nx, ny, nTriCols);
fprintf('Number of nodes    = %d\n', size(model.nodesInfo,1));
fprintf('Number of elements = %d\n', numEles);

fprintf('\nTip deflection at mid-height right node:\n');
fprintf('  FE                  = %.6e\n', uy_tip_mid);
fprintf('  Euler-Bernoulli     = %.6e\n', -w_tip_EB);
fprintf('  Timoshenko          = %.6e\n', -w_tip_Tim);
fprintf('  Rel. error vs EB    = %.6e\n', err_tip_EB_mid);
fprintf('  Rel. error vs Tim   = %.6e\n', err_tip_Tim_mid);

fprintf('\nAverage tip deflection on right edge:\n');
fprintf('  FE                  = %.6e\n', uy_tip_avg);
fprintf('  Euler-Bernoulli     = %.6e\n', -w_tip_EB);
fprintf('  Timoshenko          = %.6e\n', -w_tip_Tim);
fprintf('  Rel. error vs EB    = %.6e\n', err_tip_EB_avg);
fprintf('  Rel. error vs Tim   = %.6e\n', err_tip_Tim_avg);

fprintf('\nReaction check:\n');
fprintf('  Total vertical reaction FE = %.6e\n', Ry_total);
fprintf('  Applied load magnitude     = %.6e\n', P);
fprintf('  Relative reaction error    = %.6e\n', err_reaction);

fprintf('\nBeam-theory bending stresses at fixed end:\n');
fprintf('  sigma_top_exact = %.6e\n', sigma_top_exact);
fprintf('  sigma_bot_exact = %.6e\n', sigma_bot_exact);

%% Plots
plotMesh(model, 'Refined Cantilever Beam Mixed Mesh', false, false);

uMag = sqrt(nodal.ux.^2 + nodal.uy.^2);
plotNodalScalar(model, uMag, 'Displacement Magnitude');

plotNodalScalar(model, stressN(:,1), '\sigma_{xx}');
plotNodalScalar(model, vmN, 'Von Mises Stress');

plotDeformedMeshContour(model, result.u, 2e4, 'Refined Cantilever Deformation');