function [edgeMatrix, edgeVector, info] = Edge2_ThermoMech(edgeNodesInfo, bc, problem, loads, state)
% 2-node edge element for thermo-mechanical BC
% DoF order per node: [ux, uy, T]

edgeMatrix = zeros(6,6);
edgeVector = zeros(6,1);

% TODO: split by bc.type
% traction / flux / convection / mixed terms

info = struct();

end