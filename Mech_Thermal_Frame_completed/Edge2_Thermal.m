function [edgeMatrix, edgeVector] = Edge2_Thermal(nodeCoords, bcVals, thickness)
% 2-node thermal boundary element
%
% nodeCoords : 2 x 2 array [x y]
% bcVals     : [M, S]
% thickness  : thickness
%
% Weak form contribution:
%   int( N' * M * N ) dGamma
%   int( N' * S ) dGamma

x1 = nodeCoords(1,1); y1 = nodeCoords(1,2);
x2 = nodeCoords(2,1); y2 = nodeCoords(2,2);

L = sqrt((x2-x1)^2 + (y2-y1)^2);

M = bcVals(1);
S = bcVals(2);

% 2-point Gauss rule on line
gp = [-1/sqrt(3), 1/sqrt(3)];
w  = [1, 1];

edgeMatrix = zeros(2,2);
edgeVector = zeros(2,1);

for igp = 1:2
    s = gp(igp);

    N = [(1-s)/2, (1+s)/2];
    J = L/2;

    edgeMatrix = edgeMatrix + (N' * M * N) * J * w(igp) * thickness;
    edgeVector = edgeVector + (N' * S) * J * w(igp) * thickness;
end

end