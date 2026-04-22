function [elementMatrix, elementVector, info] = Tri3_ThermoMech(elementNodesInfo, material, problem, loads, state)
% Tri3 thermo-mechanical element
% DoF order per node: [ux, uy, T]
%
% Unknowns are solved in a linearized one-way-coupled form:
%   mechanical block includes thermal expansion load from current element temperature
%   thermal block is standard conduction
%
% Material format example:
%   [E, nu, k, alpha]

x1 = elementNodesInfo(1,2); y1 = elementNodesInfo(1,3);
x2 = elementNodesInfo(2,2); y2 = elementNodesInfo(2,3);
x3 = elementNodesInfo(3,2); y3 = elementNodesInfo(3,3);

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

A = [1 x1 y1;
     1 x2 y2;
     1 x3 y3];
detA = det(A);
area = detA / 2;

if area <= 0
    error('Tri3 element area must be positive.');
end

b1 = y2 - y3;  c1 = x3 - x2;
b2 = y3 - y1;  c2 = x1 - x3;
b3 = y1 - y2;  c3 = x2 - x1;

Bmech = 1/(2*area) * ...
    [b1  0   b2  0   b3  0;
     0   c1  0   c2  0   c3;
     c1  b1  c2  b2  c3  b3];

Bth = 1/(2*area) * [b1 b2 b3;
                    c1 c2 c3];

Kuu = thickness * area * (Bmech' * D * Bmech);
Ktt = thickness * area * (Bth' * conductivityMatrix * Bth);

couplingLoad = zeros(6,3);
if strcmpi(mechType,'planeStress')
    thermalStressVector = D * alpha * [1; 1; 0];
else
    thermalStressVector = D * alpha * [1; 1; 0];
end

NbarT = [1/3 1/3 1/3];
for iT = 1:3
    couplingLoad(:,iT) = thickness * area * (Bmech' * thermalStressVector) * NbarT(iT);
end

elementMatrix = zeros(9,9);
elementMatrix([1 2 4 5 7 8],[1 2 4 5 7 8]) = Kuu;
elementMatrix([3 6 9],[3 6 9]) = Ktt;

% Mechanical equation includes thermal expansion contribution as -Kut * T
% This is a simple readable formulation.
elementMatrix([1 2 4 5 7 8],[3 6 9]) = -couplingLoad;

elementVector = zeros(9,1);

if isfield(loads,'bodyForce')
    bodyForce = loads.bodyForce(:);
    NbarU = [1/3 0   1/3 0   1/3 0;
             0   1/3 0   1/3 0   1/3];
    elementVector([1 2 4 5 7 8]) = elementVector([1 2 4 5 7 8]) + thickness * area * NbarU' * bodyForce;
end

if isfield(loads,'heatSource')
    heatSource = loads.heatSource;
    elementVector([3 6 9]) = elementVector([3 6 9]) + thickness * area * [1/3;1/3;1/3] * heatSource;
end

if deltaTRef ~= 0
    elementVector([1 2 4 5 7 8]) = elementVector([1 2 4 5 7 8]) + couplingLoad * deltaTRef * ones(3,1);
end

info.type = 'Tri3_ThermoMech';

end