function [eleMatrix, eleVector, info] = Quad4_Thermal(eleNodesInfo, material, problem, loads, state)
% Quad4 thermal element
% DoF order per node: [T]

x = eleNodesInfo(:,2);
y = eleNodesInfo(:,3);

k = material(1);

if isfield(problem, 'thickness')
    thickness = problem.thickness;
else
    thickness = 1.0;
end

gaussPts = [-1/sqrt(3), -1/sqrt(3);
            -1/sqrt(3),  1/sqrt(3);
             1/sqrt(3), -1/sqrt(3);
             1/sqrt(3),  1/sqrt(3)];

eleMatrix = zeros(4,4);
eleVector = zeros(4,1);

for ig = 1:4
    s = gaussPts(ig,1);
    t = gaussPts(ig,2);

    N = 0.25 * [(1-s)*(1-t), (1+s)*(1-t), (1+s)*(1+t), (1-s)*(1+t)];

    dNds = 0.25 * [-(1-t),  (1-t), (1+t), -(1+t)];
    dNdt = 0.25 * [-(1-s), -(1+s), (1+s),  (1-s)];

    J = [dNds; dNdt] * [x y];
    detJ = det(J);

    if detJ <= 0
        error('Quad4 element has non-positive Jacobian.');
    end

    dNdxy = J \ [dNds; dNdt];

    B = dNdxy;
    D = k * eye(2);

    eleMatrix = eleMatrix + B' * D * B * thickness * detJ;

    if isfield(loads, 'source')
        q = loads.source;
        eleVector = eleVector + N' * q * thickness * detJ;
    end
end

info = struct();

end