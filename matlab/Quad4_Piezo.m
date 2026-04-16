function [eleMatrix, eleVector, info] = Quad4_Piezo(eleNodesInfo, material, problem, loads, state)
% Quad4 piezoelectric element
% DoF order per node: [ux, uy, phi]

% TODO:
% 1. use 2x2 Gauss integration
% 2. compute Kuu, Kpp, Kup, Kpu
% 3. assemble into 12x12 matrix

eleMatrix = zeros(12,12);
eleVector = zeros(12,1);

info = struct();

end