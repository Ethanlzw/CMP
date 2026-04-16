function [eleMatrix, eleVector, info] = Quad4_ThermoMech(eleNodesInfo, material, problem, loads, state)
% Quad4 thermo-mechanical element
% DoF order per node: [ux, uy, T]

% TODO:
% 1. use 2x2 Gauss integration
% 2. compute mechanical block Kuu
% 3. compute thermal block KTT
% 4. compute coupling blocks
% 5. assemble into 12x12 matrix

eleMatrix = zeros(12,12);
eleVector = zeros(12,1);

info = struct();

end