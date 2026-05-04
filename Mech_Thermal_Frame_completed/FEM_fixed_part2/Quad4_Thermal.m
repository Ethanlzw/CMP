function [eleMatrix, eleVector, info] = Quad4_Thermal(eleNodesInfo, material, problem, loads, state)
% Quad4 linear thermal conduction element
% DoF order per node: [T]
% Local node order:
%   1(-1,-1), 2(1,-1), 3(1,1), 4(-1,1)

x = eleNodesInfo(:,2);
y = eleNodesInfo(:,3);

if numel(material) == 1
    kxx = material(1);
    kyy = material(1);
else
    kxx = material(1);
    kyy = material(2);
end

if isfield(problem,'thickness')
    thickness = problem.thickness;
else
    thickness = 1.0;
end

Kcond = [kxx 0;
         0   kyy];

gaussPts = [-1/sqrt(3), -1/sqrt(3);
             1/sqrt(3), -1/sqrt(3);
             1/sqrt(3),  1/sqrt(3);
            -1/sqrt(3),  1/sqrt(3)];
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

    J = [dNds*x, dNds*y;
         dNdt*x, dNdt*y];
    detJ = det(J);

    if detJ <= 0
        error('Quad4_Thermal detected non-positive detJ = %.6e. Check node order or element distortion.', detJ);
    end

    gradN = J \ [dNds; dNdt];
    B = gradN;

    eleMatrix = eleMatrix + (B' * Kcond * B) * detJ * w * thickness;

    qv = 0.0;
    if nargin >= 4 && ~isempty(loads)
        if isfield(loads,'heatSource')
            qv = loads.heatSource;
        elseif isfield(loads,'source')
            qv = loads.source;
        end
    end

    eleVector = eleVector + N * qv * detJ * w * thickness;
end

info.Kcond = Kcond;
info.type = 'Quad4_Thermal';

end
