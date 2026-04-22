% =========================================================================
%  main_mixed_mech.m  –  2D plane-stress/strain FEM with mixed Tri3/Quad4
%  Cantilever beam: left edge clamped, right edge loaded downward.
% =========================================================================
clc; clear; close all;

%% ── 1. Geometry & material ───────────────────────────────────────────────
L  = 10.0;    % beam length  [m]
H  =  1.0;    % beam height  [m]
t  =  1.0;    % thickness    [m]
E  = 210e9;   % Young's modulus  [Pa]
nu = 0.30;    % Poisson's ratio  [-]
P  = 1000.0;  % total downward tip load  [N]

%% ── 2. Mesh controls ─────────────────────────────────────────────────────
nx      = 20;   % elements along x
ny      =  4;   % elements along y
nTriCols = nx;  % columns using Tri3 (0 = all Quad4, nx = all Tri3)

assert(nTriCols >= 0 && nTriCols <= nx, 'nTriCols must be in [0, nx].');

%% ── 3. Build FE model struct ─────────────────────────────────────────────
model = struct();
model.fieldNames = {'ux', 'uy'};           % DOF labels (order matters)

% Element-set registry: {setField, elementName}
model.elesSets = {
    'setTri',  'Tri3_Mech';
    'setQuad', 'Quad4_Mech'
};

model.material = struct('E', E, 'nu', nu);
model.problem  = struct('thickness', t, 'mechType', 'planeStress');
model.loads    = struct('bodyForce', [0; 0]);

%% ── 4. Mesh ──────────────────────────────────────────────────────────────
[model.nodesInfo, model.elesInfo] = buildRectMeshMixed(L, H, nx, ny, nTriCols);

fprintf('Mesh: %d nodes | %d Tri3 | %d Quad4\n', ...
    size(model.nodesInfo,1), ...
    size(model.elesInfo.setTri,1), ...
    size(model.elesInfo.setQuad,1));

%% ── 5. Boundary conditions ───────────────────────────────────────────────
tol        = 1e-12;
x          = model.nodesInfo(:,2);
leftNodes  = find(abs(x)       < tol);
rightNodes = find(abs(x - L)   < tol);

% Essential BC: clamp left edge (ux = uy = 0)
nL = numel(leftNodes);
model.bcDef.essential.nodes  = repmat(leftNodes, 2, 1);
model.bcDef.essential.fields = [repmat({'ux'}, nL, 1); repmat({'uy'}, nL, 1)];
model.bcDef.essential.values = zeros(2*nL, 1);

% Natural BC: distribute total load P equally over right-edge nodes (uy DOF)
dofInfo = getDofInfo(model);
nR = numel(rightNodes);
F_ext = zeros(size(model.nodesInfo,1) * dofInfo.nDofPerNode, 1);
uyDofs = getNodeDof(rightNodes, dofInfo.map.uy, dofInfo.nDofPerNode);
F_ext(uyDofs) = -P / nR;
model.bcDef.forceVector = F_ext;

%% ── 6. Assemble, apply BC, solve ─────────────────────────────────────────
model.bcData              = buildBCData(model);
[K, F_int]                = buildSystem(model);
[u, reaction]             = solveSystem(K, model.bcData, F_int);

%% ── 7. Post-process ──────────────────────────────────────────────────────
fields = extractNodalFields(model, u);

% Verification
uyRight = fields.uy(rightNodes);
totalRxnY = sum(reaction(getNodeDof(leftNodes, dofInfo.map.uy, dofInfo.nDofPerNode)));
fprintf('Mean tip deflection uy  = %.6e m\n',  mean(uyRight));
fprintf('Total left reaction Fy  = %.6e N\n',  totalRxnY);

% Euler-Bernoulli reference: delta = P*L^3 / (3*E*I), I = t*H^3/12
I_ref = t * H^3 / 12;
delta_ref = P * L^3 / (3 * E * I_ref);
fprintf('Euler-Bernoulli tip def = %.6e m\n', -delta_ref);

result = struct('u', u, 'fields', fields, 'reaction', reaction);

%% ── 8. Plots ─────────────────────────────────────────────────────────────
scale = 1e5;
plotMesh(model);
plotDeformedMesh(model, result, scale);
plotDeformedMeshContour(model, result, 'uy', scale);

% ── Helper: node-level DOF indices ────────────────────────────────────────
function dofs = getNodeDof(nodeIds, localDof, nDofPerNode)
    dofs = (nodeIds - 1) * nDofPerNode + localDof;
end
