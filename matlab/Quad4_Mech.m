function [eleMatrix, eleVector, info] = Quad4_Mech(eleNodesInfo, material, problem, loads, state)
% Quad4 linear mechanical element
% DoF order per node: [ux, uy]

x = eleNodesInfo(:,2);
y = eleNodesInfo(:,3);

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

D = getElasticMatrix2D(E, nu, mechType);

gaussPts = [-1/sqrt(3), -1/sqrt(3);
            -1/sqrt(3),  1/sqrt(3);
             1/sqrt(3), -1/sqrt(3);
             1/sqrt(3),  1/sqrt(3)];
weights = [1;1;1;1];

eleMatrix = zeros(8,8);
eleVector = zeros(8,1);

for ig = 1:4
    s = gaussPts(ig,1);
    t = gaussPts(ig,2);
    w = weights(ig);

    N = 0.25 * [(1-s)*(1-t), (1+s)*(1-t), (1+s)*(1+t), (1-s)*(1+t)];

    dNds = 0.25 * [-(1-t),  (1-t), (1+t), -(1+t)];
    dNdt = 0.25 * [-(1-s), -(1+s), (1+s),  (1-s)];

    J = [dNds; dNdt] * [x y];
    detJ = det(J);

    if detJ <= 0
        error('Quad4 element has non-positive Jacobian.');
    end

    dNdxy = J \ [dNds; dNdt];

    B = zeros(3,8);
    for i = 1:4
        B(:,2*i-1:2*i) = [dNdxy(1,i), 0;
                          0, dNdxy(2,i);
                          dNdxy(2,i), dNdxy(1,i)];
    end

    eleMatrix = eleMatrix + B' * D * B * thickness * detJ * w;

    if isfield(loads, 'bodyForce')
        b = loads.bodyForce(:);

        Nmat = zeros(2,8);
        for i = 1:4
            Nmat(:,2*i-1:2*i) = [N(i), 0;
                                 0,    N(i)];
        end

        eleVector = eleVector + Nmat' * b * thickness * detJ * w;
    end
end

info = struct();
info.D = D;

end