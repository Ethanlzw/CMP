function [elementMatrix, elementVector, info] = Quad4_ThermoMech(elementNodesInfo, material, problem, loads, state)
% Quad4 thermo-mechanical element
% DoF order per node: [ux, uy, T]
%
% Material format example:
%   [E, nu, k, alpha]

xCoords = elementNodesInfo(:,2);
yCoords = elementNodesInfo(:,3);

E = material(1);
nu = material(2);
k = material(3);
alpha = material(4);

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

if isfield(problem,'deltaTRef')
    deltaTRef = problem.deltaTRef;
else
    deltaTRef = 0.0;
end

D = getElasticMatrix2D(E, nu, mechType);
conductivityMatrix = k * eye(2);

if strcmpi(mechType,'planeStress')
    thermalStressVector = D * alpha * [1; 1; 0];
else
    thermalStressVector = D * alpha * [1; 1; 0];
end

gaussPoints = [-1/sqrt(3), -1/sqrt(3);
               -1/sqrt(3),  1/sqrt(3);
                1/sqrt(3), -1/sqrt(3);
                1/sqrt(3),  1/sqrt(3)];
gaussWeights = [1; 1; 1; 1];

Kuu = zeros(8,8);
Ktt = zeros(4,4);
Kut = zeros(8,4);
Fu = zeros(8,1);
Ft = zeros(4,1);

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

    Bmech = [dNdx(1) 0        dNdx(2) 0        dNdx(3) 0        dNdx(4) 0;
             0        dNdy(1) 0        dNdy(2) 0        dNdy(3) 0        dNdy(4);
             dNdy(1) dNdx(1) dNdy(2) dNdx(2) dNdy(3) dNdx(3) dNdy(4) dNdx(4)];

    Bth = dNdxy;

    NmatU = [N(1) 0    N(2) 0    N(3) 0    N(4) 0;
             0    N(1) 0    N(2) 0    N(3) 0    N(4)];

    Kuu = Kuu + Bmech' * D * Bmech * thickness * detJ * weight;
    Ktt = Ktt + Bth' * conductivityMatrix * Bth * thickness * detJ * weight;
    Kut = Kut + Bmech' * thermalStressVector * N * thickness * detJ * weight;

    if isfield(loads,'bodyForce')
        bodyForce = loads.bodyForce(:);
        Fu = Fu + NmatU' * bodyForce * thickness * detJ * weight;
    end

    if isfield(loads,'heatSource')
        heatSource = loads.heatSource;
        Ft = Ft + N' * heatSource * thickness * detJ * weight;
    end
end

elementMatrix = zeros(12,12);
elementMatrix([1 2 4 5 7 8 10 11],[1 2 4 5 7 8 10 11]) = Kuu;
elementMatrix([3 6 9 12],[3 6 9 12]) = Ktt;
elementMatrix([1 2 4 5 7 8 10 11],[3 6 9 12]) = -Kut;

elementVector = zeros(12,1);
elementVector([1 2 4 5 7 8 10 11]) = Fu;
elementVector([3 6 9 12]) = Ft;

if deltaTRef ~= 0
    elementVector([1 2 4 5 7 8 10 11]) = elementVector([1 2 4 5 7 8 10 11]) + Kut * deltaTRef * ones(4,1);
end

info.type = 'Quad4_ThermoMech';

end