function [elementMatrix, elementVector, info] = Quad4_Thermal(elementNodesInfo, material, problem, loads, state)
% Quad4 linear thermal conduction element
% DoF order per node: [T]
%
% EXAM PART:
%   Same Jacobian / Gauss integration logic as Quad4_Mech,
%   but thermal B matrix is 2x4.

xCoords = elementNodesInfo(:,2);
yCoords = elementNodesInfo(:,3);

conductivity = material(1);

if isfield(problem,'thickness')
    thickness = problem.thickness;
else
    thickness = 1.0;
end

conductivityMatrix = conductivity * eye(2);

gaussPoints = [-1/sqrt(3), -1/sqrt(3);
               -1/sqrt(3),  1/sqrt(3);
                1/sqrt(3), -1/sqrt(3);
                1/sqrt(3),  1/sqrt(3)];
gaussWeights = [1; 1; 1; 1];

elementMatrix = zeros(4,4);
elementVector = zeros(4,1);

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
    Bth = dNdxy;

    elementMatrix = elementMatrix + Bth' * conductivityMatrix * Bth * thickness * detJ * weight;

    if isfield(loads,'heatSource')
        heatSource = loads.heatSource;
        elementVector = elementVector + N' * heatSource * thickness * detJ * weight;
    end
end

info.type = 'Quad4_Thermal';

end