function [edgeMatrix, edgeVector, info] = Edge2_Thermal(edgeNodesInfo, bc, problem, loads, state)
% 2-node edge element for thermal BC
% DoF order per node: [T]

x1 = edgeNodesInfo(1,2);  y1 = edgeNodesInfo(1,3);
x2 = edgeNodesInfo(2,2);  y2 = edgeNodesInfo(2,3);

Ledge = hypot(x2 - x1, y2 - y1);
Jedge = Ledge / 2;

gaussPts = [-1/sqrt(3); 1/sqrt(3)];
weights  = [1; 1];

edgeMatrix = zeros(2,2);
edgeVector = zeros(2,1);

for ig = 1:2
    s = gaussPts(ig);
    w = weights(ig);

    N = [(1 - s)/2, (1 + s)/2];

    switch lower(bc.type)
        case 'flux'
            qn = bc.value(1);
            edgeVector = edgeVector + N' * qn * Jedge * w;

        case 'convection'
            h = bc.value(1);
            Tinf = bc.value(2);
            edgeMatrix = edgeMatrix + (N' * N) * h * Jedge * w;
            edgeVector = edgeVector + N' * h * Tinf * Jedge * w;

        otherwise
            error('Unsupported thermal edge BC type.');
    end
end

info = struct();
info.length = Ledge;

end