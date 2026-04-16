function [edgeMatrix, edgeVector, info] = Edge2_MechTraction(edgeNodesInfo, bc, problem, loads, state)
% 2-node edge element for mechanical traction
% DoF order per node: [ux, uy]

x1 = edgeNodesInfo(1,2);  y1 = edgeNodesInfo(1,3);
x2 = edgeNodesInfo(2,2);  y2 = edgeNodesInfo(2,3);

Ledge = hypot(x2 - x1, y2 - y1);
Jedge = Ledge / 2;

gaussPts = [-1/sqrt(3); 1/sqrt(3)];
weights  = [1; 1];

edgeMatrix = zeros(4,4);
edgeVector = zeros(4,1);

for ig = 1:2
    s = gaussPts(ig);
    w = weights(ig);

    N1 = (1 - s) / 2;
    N2 = (1 + s) / 2;

    N = [N1 0 N2 0;
         0 N1 0 N2];

    traction = bc.value;
    edgeVector = edgeVector + N' * traction * Jedge * w;
end

info = struct();
info.length = Ledge;

end