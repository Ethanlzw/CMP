function [eleMatrix, eleVector, info] = Tri3_ThermoMech(eleNodesInfo, material, problem, loads, state)
% Tri3 linear thermomechanical element
% DoF order per node: [ux, uy, T]
%
% Element block system:
%   [Kuu KuT] [u] = [fu]
%   [ 0  KTT] [T]   [fT]
%
% The coupling is one-way: temperature creates thermal strain, but the
% mechanical displacement does not feed back into the steady thermal field.

x1 = eleNodesInfo(1,2); y1 = eleNodesInfo(1,3);
x2 = eleNodesInfo(2,2); y2 = eleNodesInfo(2,3);
x3 = eleNodesInfo(3,2); y3 = eleNodesInfo(3,3);

E     = material(1);
nu    = material(2);
alpha = material(3);
kxx   = material(4);
kyy   = material(5);

if isfield(problem,'thickness')
    thickness = problem.thickness;
else
    thickness = 1.0;
end

if isfield(problem,'mechType')
    mechType = problem.mechType;
else
    mechType = 'planeStress';
end

if isfield(problem,'Tref')
    Tref = problem.Tref;
else
    Tref = 0.0;
end

D = getElasticMatrix2D(E, nu, mechType);
Kcond = [kxx 0;
         0   kyy];

A = [1 x1 y1;
     1 x2 y2;
     1 x3 y3];
detA = det(A);
Area = detA / 2;

if Area <= 0
    error('Tri3_ThermoMech has non-positive area. Check node ordering.');
end

b1 = y2 - y3;  b2 = y3 - y1;  b3 = y1 - y2;
c1 = x3 - x2;  c2 = x1 - x3;  c3 = x2 - x1;

Bu = [b1 0  b2 0  b3 0;
      0  c1 0  c2 0  c3;
      c1 b1 c2 b2 c3 b3] / detA;

BT = [b1 b2 b3;
      c1 c2 c3] / detA;

Kuu = thickness * Area * (Bu' * D * Bu);
KTT = thickness * Area * (BT' * Kcond * BT);

epsTh = alpha * [1; 1; 0];
Nint = Area * [1/3, 1/3, 1/3];
KuT = -thickness * (Bu' * D * epsTh) * Nint;

fMech = zeros(6,1);
fTherm = zeros(3,1);

% Reference temperature contribution. This is needed when Tref is not zero.
fMech = fMech - thickness * Area * (Bu' * D * epsTh) * Tref;

if nargin >= 4 && ~isempty(loads)
    if isfield(loads,'bodyForce')
        bodyForce = loads.bodyForce(:);
        if numel(bodyForce) ~= 2
            error('Thermomechanical bodyForce must be [bx; by].');
        end
        bx = bodyForce(1);
        by = bodyForce(2);
        fMech = fMech + thickness * Area / 3 * [bx; by; bx; by; bx; by];
    end

    if isfield(loads,'heatSource')
        qv = loads.heatSource;
        fTherm = fTherm + thickness * Area / 3 * qv * [1; 1; 1];
    end
end

mechDofs = [1 2 4 5 7 8];
tempDofs = [3 6 9];

eleMatrix = zeros(9,9);
eleMatrix(mechDofs, mechDofs) = Kuu;
eleMatrix(mechDofs, tempDofs) = KuT;
eleMatrix(tempDofs, tempDofs) = KTT;

eleVector = zeros(9,1);
eleVector(mechDofs) = fMech;
eleVector(tempDofs) = fTherm;

info.Kuu = Kuu;
info.KuT = KuT;
info.KTT = KTT;
info.Bu = Bu;
info.BT = BT;
info.Area = Area;
info.type = 'Tri3_ThermoMech';

end
