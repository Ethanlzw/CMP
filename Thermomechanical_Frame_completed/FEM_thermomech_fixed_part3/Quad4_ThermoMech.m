function [eleMatrix, eleVector, info] = Quad4_ThermoMech(eleNodesInfo, material, problem, loads, state)
% Quad4 linear thermomechanical element
% DoF order per node: [ux, uy, T]
%
% The local node order is assumed to be:
%   1: bottom-left, 2: bottom-right, 3: top-right, 4: top-left
% with natural coordinates:
%   1: (-1,-1), 2: (1,-1), 3: (1,1), 4: (-1,1)

x = eleNodesInfo(:,2);
y = eleNodesInfo(:,3);

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

gaussPts = [-1/sqrt(3), -1/sqrt(3);
             1/sqrt(3), -1/sqrt(3);
             1/sqrt(3),  1/sqrt(3);
            -1/sqrt(3),  1/sqrt(3)];
weights = [1; 1; 1; 1];

Kuu = zeros(8,8);
KuT = zeros(8,4);
KTT = zeros(4,4);

fMech = zeros(8,1);
fTherm = zeros(4,1);

epsTh = alpha * [1; 1; 0];

for igp = 1:4
    s = gaussPts(igp,1);
    t = gaussPts(igp,2);
    w = weights(igp);

    N = 0.25 * [(1-s)*(1-t);
                (1+s)*(1-t);
                (1+s)*(1+t);
                (1-s)*(1+t)];

    dNds = 0.25 * [-(1-t);
                    (1-t);
                    (1+t);
                   -(1+t)];

    dNdt = 0.25 * [-(1-s);
                   -(1+s);
                    (1+s);
                    (1-s)];

    J = [dNds'*x, dNds'*y;
         dNdt'*x, dNdt'*y];

    detJ = det(J);
    if detJ <= 0
        error('Quad4_ThermoMech has non-positive detJ. Check node ordering or element distortion.');
    end

    gradN = J \ [dNds'; dNdt'];
    dNdx = gradN(1,:)';
    dNdy = gradN(2,:)';

    Bu = [dNdx(1) 0       dNdx(2) 0       dNdx(3) 0       dNdx(4) 0;
          0       dNdy(1) 0       dNdy(2) 0       dNdy(3) 0       dNdy(4);
          dNdy(1) dNdx(1) dNdy(2) dNdx(2) dNdy(3) dNdx(3) dNdy(4) dNdx(4)];

    BT = [dNdx';
          dNdy'];

    Kuu = Kuu + (Bu' * D * Bu) * detJ * w * thickness;
    KTT = KTT + (BT' * Kcond * BT) * detJ * w * thickness;
    KuT = KuT - (Bu' * D * epsTh) * N' * detJ * w * thickness;

    % Reference temperature contribution. This is needed when Tref is not zero.
    fMech = fMech - (Bu' * D * epsTh) * Tref * detJ * w * thickness;

    if nargin >= 4 && ~isempty(loads)
        if isfield(loads,'bodyForce')
            bodyForce = loads.bodyForce(:);
            if numel(bodyForce) ~= 2
                error('Thermomechanical bodyForce must be [bx; by].');
            end
            Nmat = [N(1) 0    N(2) 0    N(3) 0    N(4) 0;
                    0    N(1) 0    N(2) 0    N(3) 0    N(4)];
            fMech = fMech + Nmat' * bodyForce * detJ * w * thickness;
        end

        if isfield(loads,'heatSource')
            qv = loads.heatSource;
            fTherm = fTherm + N * qv * detJ * w * thickness;
        end
    end
end

eleMatrix = zeros(12,12);

mechDofs = [1 2 4 5 7 8 10 11];
tempDofs = [3 6 9 12];

eleMatrix(mechDofs, mechDofs) = Kuu;
eleMatrix(mechDofs, tempDofs) = KuT;
eleMatrix(tempDofs, tempDofs) = KTT;

eleVector = zeros(12,1);
eleVector(mechDofs) = fMech;
eleVector(tempDofs) = fTherm;

info.Kuu = Kuu;
info.KuT = KuT;
info.KTT = KTT;
info.type = 'Quad4_ThermoMech';

end
