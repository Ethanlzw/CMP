function [elementMatrix, elementVector, info] = Tri3_Thermal(elementNodesInfo, material, problem, loads, state)
% Tri3 linear thermal conduction element
% DoF order per node: [T]
%
% Governing form:
%   K = integral( Bth' * k * Bth * thickness * dA )

x1 = elementNodesInfo(1,2); y1 = elementNodesInfo(1,3);
x2 = elementNodesInfo(2,2); y2 = elementNodesInfo(2,3);
x3 = elementNodesInfo(3,2); y3 = elementNodesInfo(3,3);

conductivity = material(1);

if isfield(problem,'thickness')
    thickness = problem.thickness;
else
    thickness = 1.0;
end

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

Bth = 1/(2*area) * [b1 b2 b3;
                    c1 c2 c3];

conductivityMatrix = conductivity * eye(2);

elementMatrix = thickness * area * (Bth' * conductivityMatrix * Bth);

elementVector = zeros(3,1);
if isfield(loads,'heatSource')
    heatSource = loads.heatSource;
    Nbar = [1/3; 1/3; 1/3];
    elementVector = elementVector + thickness * area * Nbar * heatSource;
end

info.area = area;
info.Bth = Bth;
info.type = 'Tri3_Thermal';

end