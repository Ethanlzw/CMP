function [elementMatrix, elementVector, info] = Tri3_Mech(elementNodesInfo, material, problem, loads, state)
% Tri3 linear mechanical element
% DoF order per node: [ux, uy]
%
% EXAM PART:
%   This file contains high-frequency exam content:
%   - area of triangle
%   - B matrix of Tri3
%   - K = t*A*B'*D*B

x1 = elementNodesInfo(1,2); y1 = elementNodesInfo(1,3);
x2 = elementNodesInfo(2,2); y2 = elementNodesInfo(2,3);
x3 = elementNodesInfo(3,2); y3 = elementNodesInfo(3,3);

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

B = 1/(2*area) * ...
    [b1  0   b2  0   b3  0;
     0   c1  0   c2  0   c3;
     c1  b1  c2  b2  c3  b3];

elementMatrix = thickness * area * (B' * D * B);

elementVector = zeros(6,1);
if isfield(loads,'bodyForce')
    bodyForce = loads.bodyForce(:);
    Nbar = [1/3 0   1/3 0   1/3 0;
            0   1/3 0   1/3 0   1/3];
    elementVector = elementVector + thickness * area * Nbar' * bodyForce;
end

info.area = area;
info.B = B;
info.D = D;
info.type = 'Tri3_Mech';

end