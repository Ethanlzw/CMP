function [Ke, Fe, info] = Tri3_Mech(eleNodesInfo, material, problem, loads, ~)
% Tri3_Mech  Constant-strain triangle (CST) element for 2D mechanics.
%
%   DOF order per node: [ux, uy]
%   Uses analytic B-matrix; element is constant-strain (one Gauss point).
%
%   Inputs:
%     eleNodesInfo – [3×3] columns: [nodeId, x, y]
%     material     – struct with fields E, nu
%     problem      – struct with fields thickness, mechType
%     loads        – struct with optional field bodyForce [2×1]
%
%   Outputs:
%     Ke   – [6×6] element stiffness matrix
%     Fe   – [6×1] element load vector
%     info – struct with diagnostic fields (area, B, D)

x = eleNodesInfo(:,2);
y = eleNodesInfo(:,3);

E  = material.E;
nu = material.nu;
t  = getfield_default(problem, 'thickness', 1.0);
mt = getfield_default(problem, 'mechType',  'planeStress');

% ── Element area ──────────────────────────────────────────────────────────
A2 = det([1 x(1) y(1);
          1 x(2) y(2);
          1 x(3) y(3)]);
A = A2 / 2;
assert(A > 0, 'Tri3_Mech: element has non-positive area (check node ordering).');

% ── Strain-displacement matrix B (constant over element) ──────────────────
%   Derived from shape functions N_i = (a_i + b_i*x + c_i*y) / (2A)
b = [y(2)-y(3),  y(3)-y(1),  y(1)-y(2)];   % dN/dx * 2A
c = [x(3)-x(2),  x(1)-x(3),  x(2)-x(1)];   % dN/dy * 2A

B = (1/(2*A)) * ...
    [b(1), 0,    b(2), 0,    b(3), 0;
     0,    c(1), 0,    c(2), 0,    c(3);
     c(1), b(1), c(2), b(2), c(3), b(3)];

% ── Constitutive matrix ───────────────────────────────────────────────────
D  = getElasticMatrix2D(E, nu, mt);
Ke = (t * A) * (B' * D * B);

% ── Load vector (body force, centroid quadrature) ─────────────────────────
Fe = zeros(6, 1);
if isfield(loads, 'bodyForce')
    Nbar = (1/3) * [1,0, 1,0, 1,0;
                    0,1, 0,1, 0,1];   % shape functions at centroid
    Fe = (t * A) * Nbar' * loads.bodyForce(:);
end

info = struct('area', A, 'B', B, 'D', D);

end

% ── Local helper ─────────────────────────────────────────────────────────
function v = getfield_default(s, fname, default)
    if isfield(s, fname), v = s.(fname); else, v = default; end
end
