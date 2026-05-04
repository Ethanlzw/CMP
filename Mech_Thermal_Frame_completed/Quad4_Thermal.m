function [eleMatrix, eleVector, info] = Quad4_Thermal(eleNodesInfo, material, problem, loads, state)
% Quad4 linear thermal conduction element
% DoF order per node: [T]

x = eleNodesInfo(:,2);
y = eleNodesInfo(:,3);

kxx = material(1);
kyy = material(2);

if isfield(problem,'thickness')
    thickness = problem.thickness;
else
    thickness = 1.0;
end

Kcond = [kxx 0;
         0   kyy];

gaussPts = [-1/sqrt(3), -1/sqrt(3);
            -1/sqrt(3),  1/sqrt(3);
             1/sqrt(3), -1/sqrt(3);
             1/sqrt(3),  1/sqrt(3)];
weights = [1; 1; 1; 1];

eleMatrix = zeros(4,4);
eleVector = zeros(4,1);

for igp = 1:4
    s = gaussPts(igp,1);
    t = gaussPts(igp,2);
    w = weights(igp);
    
    N = 0.25 * [(1-s)*(1-t);
                (1+s)*(1-t);
                (1+s)*(1+t);
                (1-s)*(1+t)];

    dNds = 0.25 * [-(1-t),  (1-t),  (1+t), -(1+t)];
    dNdt = 0.25 * [-(1-s), -(1+s),  (1+s),  (1-s)];
%% Complete the following section
    J = [dNds'*x, dNds'*y;
         dNdt'*x, dNdt'*y];

    detJ = det(J);
    invJ = inv(J);

    gradN = invJ * [dNds'; dNdt'];   % 2 x 4
    B = gradN;

    eleMatrix = eleMatrix + (B' * Kcond * B) * detJ * w * thickness;

    if nargin >= 4 && ~isempty(loads)
        if isfield(loads,'heatSource')
            qv = loads.heatSource;
            eleVector = eleVector + N * qv * detJ * w * thickness;
        end
    end
end

info.Kcond = Kcond;
info.type = 'Quad4_Thermal';

end