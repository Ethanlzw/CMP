function [eleMatrix, eleVector, info] = Tri3_Thermal(eleNodesInfo, material, problem, loads, state)
% Tri3 linear thermal conduction element
% DoF order per node: [T]

x1 = eleNodesInfo(1,2); y1 = eleNodesInfo(1,3);
x2 = eleNodesInfo(2,2); y2 = eleNodesInfo(2,3);
x3 = eleNodesInfo(3,2); y3 = eleNodesInfo(3,3);

kxx = material(1);
kyy = material(2);

if isfield(problem,'thickness')
    thickness = problem.thickness;
else
    thickness = 1.0;
end

A = [1 x1 y1;
     1 x2 y2;
     1 x3 y3];
detA = det(A);
Area = detA / 2;

b1 = y2 - y3;  b2 = y3 - y1;  b3 = y1 - y2;
c1 = x3 - x2;  c2 = x1 - x3;  c3 = x2 - x1;

%% Complete the following section
B = [b1 b2 b3;
     c1 c2 c3] / detA;

Kcond = [kxx 0;
         0   kyy];

eleMatrix = thickness * Area * (B' * Kcond * B);

eleVector = zeros(3,1);

if nargin >= 4 && ~isempty(loads)
    if isfield(loads,'heatSource')
        qv = loads.heatSource;
        eleVector = eleVector + thickness * Area * [1/3; 1/3; 1/3] * qv; % Complete this line
    end
end
%%
info.Area = Area;
info.B = B;
info.Kcond = Kcond;
info.type = 'Tri3_Thermal';

end