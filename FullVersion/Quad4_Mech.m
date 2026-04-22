function [elementMatrix, elementVector, info] = Quad4_Mech(elementNodesInfo, material, problem, loads, state)
% Quad4 linear mechanical element
% DoF order per node: [ux, uy]
%
% EXAM PART:
%   This file contains high-frequency exam content:
%   - shape function derivatives in (s,t)
%   - Jacobian matrix
%   - B matrix of Quad4
%   - 2x2 Gauss integration

xCoords = elementNodesInfo(:,2);
yCoords = elementNodesInfo(:,3);

E = material(1);
nu = material(2);

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

D = getElasticMatrix2D(E, nu, mechType);

gaussPoints = [-1/sqrt(3), -1/sqrt(3);
               -1/sqrt(3),  1/sqrt(3);
                1/sqrt(3), -1/sqrt(3);
                1/sqrt(3),  1/sqrt(3)];
gaussWeights = [1; 1; 1; 1];

elementMatrix = zeros(8,8);
elementVector = zeros(8,1);

for iGauss = 1:4
    s = gaussPoints(iGauss,1);
    t = gaussPoints(iGauss,2);
    weight = gaussWeights(iGauss);

    N = 0.25 * [(1-s)*(1-t), (1+s)*(1-t), (1+s)*(1+t), (1-s)*(1+t)];

    dNds = 0.25 * [-(1-t),  (1-t), (1+t), -(1+t)];
    dNdt = 0.25 * [-(1-s), -(1+s), (1+s),  (1-s)];

    J = [dNds; dNdt] * [xCoords yCoords];
    detJ = det(J);

    if detJ <= 0
        error('Quad4 element has non-positive Jacobian.');
    end

    dNdxy = J \ [dNds; dNdt];
    dNdx = dNdxy(1,:);
    dNdy = dNdxy(2,:);

    B = [dNdx(1) 0        dNdx(2) 0        dNdx(3) 0        dNdx(4) 0;
         0        dNdy(1) 0        dNdy(2) 0        dNdy(3) 0        dNdy(4);
         dNdy(1) dNdx(1) dNdy(2) dNdx(2) dNdy(3) dNdx(3) dNdy(4) dNdx(4)];

    elementMatrix = elementMatrix + B' * D * B * thickness * detJ * weight;

    if isfield(loads,'bodyForce')
        bodyForce = loads.bodyForce(:);
        Nmat = [N(1) 0    N(2) 0    N(3) 0    N(4) 0;
                0    N(1) 0    N(2) 0    N(3) 0    N(4)];
        elementVector = elementVector + Nmat' * bodyForce * thickness * detJ * weight;
    end
end

info.D = D;
info.type = 'Quad4_Mech';

end