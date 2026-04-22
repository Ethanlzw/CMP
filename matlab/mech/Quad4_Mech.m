function [Ke, Fe, info] = Quad4_Mech(eleNodesInfo, material, problem, loads, ~)
% Quad4_Mech  Bilinear quadrilateral element for 2D mechanics.
%
%   DOF order per node: [ux, uy]
%   Uses 2×2 Gauss quadrature in (s,t) ∈ [-1,1]².
%
%   Inputs:
%     eleNodesInfo – [4×3] columns: [nodeId, x, y]  (counter-clockwise)
%     material     – struct with fields E, nu
%     problem      – struct with fields thickness, mechType
%     loads        – struct with optional field bodyForce [2×1]
%
%   Outputs:
%     Ke   – [8×8] element stiffness matrix
%     Fe   – [8×1] element load vector
%     info – struct with diagnostic field D

x = eleNodesInfo(:,2);
y = eleNodesInfo(:,3);

E  = material.E;
nu = material.nu;
t  = getfield_default(problem, 'thickness', 1.0);
mt = getfield_default(problem, 'mechType',  'planeStress');

D  = getElasticMatrix2D(E, nu, mt);
Ke = zeros(8, 8);
Fe = zeros(8, 1);

% ── 2×2 Gauss points and weights ─────────────────────────────────────────
gp = 1/sqrt(3);
gaussPts = [-gp, -gp;  gp, -gp;  gp, gp;  -gp, gp];
w        = ones(4, 1);

for ig = 1:4
    s = gaussPts(ig,1);
    tt = gaussPts(ig,2);   % 't' is reserved in MATLAB; use 'tt'

    % Shape functions and their natural-coord derivatives
    N    = 0.25 * [(1-s)*(1-tt), (1+s)*(1-tt), (1+s)*(1+tt), (1-s)*(1+tt)];
    dNds = 0.25 * [-(1-tt),  (1-tt),  (1+tt), -(1+tt)];
    dNdt = 0.25 * [-(1-s),  -(1+s),   (1+s),   (1-s) ];

    % Jacobian and its inverse
    J    = [dNds; dNdt] * [x, y];
    detJ = det(J);
    assert(detJ > 0, 'Quad4_Mech: non-positive Jacobian (check node ordering).');

    dNdxy = J \ [dNds; dNdt];  % [2×4] physical derivatives

    % Strain-displacement matrix B [3×8]
    B = zeros(3, 8);
    for i = 1:4
        B(:, 2*i-1:2*i) = [dNdxy(1,i), 0;
                            0,          dNdxy(2,i);
                            dNdxy(2,i), dNdxy(1,i)];
    end

    fac = t * detJ * w(ig);
    Ke  = Ke + fac * (B' * D * B);

    if isfield(loads, 'bodyForce')
        Nmat = zeros(2, 8);
        for i = 1:4
            Nmat(:, 2*i-1:2*i) = [N(i), 0; 0, N(i)];
        end
        Fe = Fe + fac * (Nmat' * loads.bodyForce(:));
    end
end

info = struct('D', D);

end

% ── Local helper ─────────────────────────────────────────────────────────
function v = getfield_default(s, fname, default)
    if isfield(s, fname), v = s.(fname); else, v = default; end
end
