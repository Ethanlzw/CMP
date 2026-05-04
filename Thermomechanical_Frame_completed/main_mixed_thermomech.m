clc;
clear;
close all;

%% Geometry and material
L = 2.0;
H = 1.0;

E = 210e9;
nu = 0.30;
alpha = 1.2e-5;
kxx = 10.0;
kyy = 10.0;
thickness = 1.0;

Tref = 0.0;
DeltaT = 100.0;
Tuni = Tref + DeltaT;

%% Mesh controls
nx = 16;
ny = 8;
nTriCols = nx;  

if nTriCols < 0 || nTriCols > nx
    error('nTriCols must satisfy 0 <= nTriCols <= nx.');
end

%% Create FE model
model = struct();
model.fieldNames = {'ux','uy','T'};

model.elesType = {
    1, 'Tri3_ThermoMech',  3, 3;
    2, 'Quad4_ThermoMech', 3, 4
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

node = @(i,j) j*(nx+1) + i + 1;

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
model.problem.mechType = 'planeStress';
model.problem.Tref = Tref;

%% Section/material data
numEles = eleID;
model.section = zeros(numEles, 6);   % [eleID, E, nu, alpha, kxx, kyy]
model.section(:,1) = 1:numEles;
for e = 1:numEles
    model.section(e,2:6) = [E, nu, alpha, kxx, kyy];
end

%% Element loads
model.eleLoadData = cell(numEles,1);
for e = 1:numEles
    model.eleLoadData{e} = struct('bodyForce', [0;0], 'heatSource', 0.0);
end

%% Global nodal vector
nNodeDoF = model.elesType{1,3};
model.nodalVector = zeros(size(model.nodesInfo,1)*nNodeDoF, 1);

%% Essential BCs
% Uniform temperature on all boundary nodes
% Minimal displacement constraints:
%   node (0,0): ux=0, uy=0
%   node (L,0): uy=0
%
% Since buildBCData applies all DoFs at a node, this framework is currently
% node-group based, not component-wise. So for this validation, use:
%   left-bottom node: [0,0,Tuni]
%   right-bottom node: [alpha*DeltaT*L, 0, Tuni]
% and prescribe Tuni on all other boundary nodes consistently.
%
% This still matches the exact solution and avoids rigid body motion.

essBC = {};

% All boundary nodes get exact displacement+temperature if on constrained corners,
% otherwise exact temperature only is not currently supported component-wise.
% So for this benchmark, prescribe the exact full solution on all boundary nodes.
%
% This is effectively a thermomechanical patch test with exact boundary field.

for j = 0:ny
    for i = 0:nx
        if i==0 || i==nx || j==0 || j==ny
            nid = node(i,j);
            x = model.nodesInfo(nid,2);
            y = model.nodesInfo(nid,3);

            ux = alpha * DeltaT * x;
            uy = alpha * DeltaT * y;

            essBC(end+1,1:2) = {nid, [ux, uy, Tuni]};    %#ok<SAGROW>
        end
    end
end

model.essBC = essBC;

% No natural BCs needed
model.natBC = {};

%% Solve FE problem
[globalMatrix, globalVector] = buildSystem(model);
bcData = buildBCData(model);
[reaction, solution] = solveSystem(globalMatrix, globalVector, bcData);

result.u = solution;
result.f = reaction;

%% Postprocessing
nodal = extractNodalFields(model, result.u);
results = postprocessModel(model, result.u);

%% Exact nodal fields
xN = model.nodesInfo(:,2);
yN = model.nodesInfo(:,3);

ux_exact = alpha * DeltaT * xN;
uy_exact = alpha * DeltaT * yN;
T_exact  = Tuni * ones(size(xN));

ux_err = nodal.ux - ux_exact;
uy_err = nodal.uy - uy_exact;
T_err  = nodal.T  - T_exact;

u_err_max = max([max(abs(ux_err)), max(abs(uy_err))]);
u_err_L2  = sqrt(mean([ux_err; uy_err].^2));
T_err_max = max(abs(T_err));
T_err_L2  = sqrt(mean(T_err.^2));

%% Element comparisons
numEles = numel(results);

strain_FE = zeros(numEles,3);
strainTh_FE = zeros(numEles,3);
stress_FE = zeros(numEles,3);

strain_exact = repmat([alpha*DeltaT, alpha*DeltaT, 0], numEles, 1);
strainTh_exact = strain_exact;
stress_exact = zeros(numEles,3);

for e = 1:numEles
    strain_FE(e,:) = results(e).value.strain(:)';
    strainTh_FE(e,:) = results(e).value.thermalStrain(:)';
    stress_FE(e,:) = results(e).value.stress(:)';
end

strain_err_max = max(abs(strain_FE - strain_exact), [], 'all');
strainTh_err_max = max(abs(strainTh_FE - strainTh_exact), [], 'all');
stress_err_max = max(abs(stress_FE - stress_exact), [], 'all');

%% Nodal averaged element fields
vmN = averageFieldToNodes(model, results, 'vonMises');
stressN = averageFieldToNodes(model, results, 'stress');

%% Print results
disp('===== Thermomechanical Free Expansion Validation =====');
fprintf('nx = %d, ny = %d, nTriCols = %d\n', nx, ny, nTriCols);
fprintf('Number of nodes    = %d\n', size(model.nodesInfo,1));
fprintf('Number of elements = %d\n', numEles);

fprintf('\nNodal errors:\n');
fprintf('  Max displacement error = %.6e\n', u_err_max);
fprintf('  L2 displacement error  = %.6e\n', u_err_L2);
fprintf('  Max temperature error  = %.6e\n', T_err_max);
fprintf('  L2 temperature error   = %.6e\n', T_err_L2);

fprintf('\nElement errors:\n');
fprintf('  Max strain error         = %.6e\n', strain_err_max);
fprintf('  Max thermal strain error = %.6e\n', strainTh_err_max);
fprintf('  Max stress error         = %.6e\n', stress_err_max);

%% Plot
plotMesh(model, 'Thermomechanical Free Expansion Mesh', false, false);

plotNodalScalar(model, nodal.T, 'Temperature');
plotNodalScalar(model, nodal.ux, 'u_x');
plotNodalScalar(model, nodal.uy, 'u_y');
plotNodalScalar(model, vmN, 'Von Mises Stress');

plotNodalScalar(model, stressN(:,1), '\sigma_{xx}');
plotNodalScalar(model, stressN(:,2), '\sigma_{yy}');
plotNodalScalar(model, stressN(:,3), '\tau_{xy}');

plotDeformedMeshContour(model, result.u, 50, 'Free Thermal Expansion');