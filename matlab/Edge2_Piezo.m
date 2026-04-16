function [edgeMatrix, edgeVector, info] = Edge2_Piezo(edgeNodesInfo, bc, problem, loads, state)
% 2-node edge element for piezoelectric BC
% DoF order per node: [ux, uy, phi]

edgeMatrix = zeros(6,6);
edgeVector = zeros(6,1);

% TODO: split by bc.type
% traction / surface charge / electric flux

info = struct();

end