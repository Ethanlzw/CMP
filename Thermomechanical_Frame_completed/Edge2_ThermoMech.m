function [edgeMatrix, edgeVector] = Edge2_ThermoMech(nodeCoords, bcVals, thickness)
% 2-node boundary element for thermomechanical problems
%
% nodeCoords : 2 x 2 array [x y]
% bcVals     : [tx, ty, M, S]
% thickness  : thickness
%
% nodal DoF order on edge:
%   [ux1 uy1 T1 ux2 uy2 T2]'

x1 = nodeCoords(1,1); y1 = nodeCoords(1,2);
x2 = nodeCoords(2,1); y2 = nodeCoords(2,2);

L = sqrt((x2-x1)^2 + (y2-y1)^2);

tx = bcVals(1);
ty = bcVals(2);
M  = bcVals(3);
S  = bcVals(4);

gp = [-1/sqrt(3), 1/sqrt(3)];
w  = [1, 1];

edgeMatrix = zeros(6,6);
edgeVector = zeros(6,1);

for igp = 1:2
    s = gp(igp);

    N = [(1-s)/2, (1+s)/2];
    J = L/2;


    %% Complete the following session

    % Mechanical traction contribution
    Nmat_u = [N(1) 0    N(2) 0;
              0    N(1) 0    N(2)];
    fu = Nmat_u' * [tx; ty] * J * w(igp) * thickness;

    % Thermal Robin / source contribution
    KT = (N' * M * N) * J * w(igp) * thickness;
    fT = (N' * S) * J * w(igp) * thickness;

    % Assemble into [ux1 uy1 T1 ux2 uy2 T2]
    mechDofs = [1 2 4 5];
    tempDofs = [3 6];

    edgeVector(mechDofs) = edgeVector(mechDofs) + fu;
    edgeVector(tempDofs) = edgeVector(tempDofs) + fT;
    edgeMatrix(tempDofs, tempDofs) = edgeMatrix(tempDofs, tempDofs) + KT;
end

end