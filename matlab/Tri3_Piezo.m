function [eleMatrix, eleVector, info] = Tri3_Piezo(eleNodesInfo, material, problem, loads, state)
% Tri3 piezoelectric element
% DoF order per node: [ux, uy, phi]

% TODO:
% 1. build mechanical stiffness Kuu
% 2. build dielectric stiffness Kpp
% 3. build piezoelectric coupling Kup / Kpu
% 4. assemble into 9x9 matrix

eleMatrix = zeros(9,9);
eleVector = zeros(9,1);

info = struct();

end