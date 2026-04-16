function [eleMatrix, eleVector, info] = Tri3_Mech(eleNodesInfo, material, problem, loads, state)
% Tri3 linear mechanical element
% DoF order per node: [ux, uy]

x1 = eleNodesInfo(1,2); y1 = eleNodesInfo(1,3);
x2 = eleNodesInfo(2,2); y2 = eleNodesInfo(2,3);
x3 = eleNodesInfo(3,2); y3 = eleNodesInfo(3,3);

E  = material(1);
nu = material(2);

if isfield(problem, 'thickness')
    thickness = problem.thickness;
else
    thickness = 1.0;
end

if isfield(problem, 'mechType')
    mechType = problem.mechType;
else
    mechType = 'planeStress';
end

A2 = det([1 x1 y1;
          1 x2 y2;
          1 x3 y3]);

A = A2 / 2;

if A <= 0
    error('Tri3 element area must be positive.');
end

b1 = y2 - y3;  c1 = x3 - x2;
b2 = y3 - y1;  c2 = x1 - x3;
b3 = y1 - y2;  c3 = x2 - x1;

B = 1 / (2*A) * ...
    [b1  0  b2  0  b3  0;
     0  c1  0  c2  0  c3;
     c1 b1 c2 b2 c3 b3];

D = getElasticMatrix2D(E, nu, mechType);

eleMatrix = thickness * A * (B' * D * B);
eleVector = zeros(6,1);

if isfield(loads, 'bodyForce')
    b = loads.bodyForce(:);
    Nbar = [1/3 0 1/3 0 1/3 0;
            0 1/3 0 1/3 0 1/3];
    eleVector = eleVector + thickness * A * Nbar' * b;
end

info = struct();
info.area = A;
info.B = B;
info.D = D;

end