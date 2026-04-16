function [eleMatrix, eleVector, info] = Tri3_ThermoMech(eleNodesInfo, material, problem, loads, state)
% Tri3 thermo-mechanical element
% DoF order per node: [ux, uy, T]

% TODO:
% 1. build Kuu from mechanics
% 2. build KTT from thermal
% 3. build coupling blocks KuT / KTu if needed
% 4. assemble into 9x9 block matrix

eleMatrix = zeros(9,9);
eleVector = zeros(9,1);

info = struct();

end