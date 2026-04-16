function [eleMatrix, eleVector, info] = Tri3_Thermal(eleNodesInfo, material, problem, loads, state)
% Tri3 thermal element
% DoF order per node: [T]

x1 = eleNodesInfo(1,2); y1 = eleNodesInfo(1,3);
x2 = eleNodesInfo(2,2); y2 = eleNodesInfo(2,3);
x3 = eleNodesInfo(3,2); y3 = eleNodesInfo(3,3);

k = material(1);

if isfield(problem, 'thickness')
    thickness = problem.thickness;
else
    thickness = 1.0;
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

B = 1 / (2*A) * [b1 b2 b3;
                 c1 c2 c3];

D = k * eye(2);

eleMatrix = thickness * A * (B' * D * B);
eleVector = zeros(3,1);

if isfield(loads, 'source')
    q = loads.source;
    eleVector = eleVector + thickness * A * q / 3 * ones(3,1);
end

info = struct();
info.area = A;

end