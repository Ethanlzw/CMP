function edgeVector = Edge2_MechTraction(nodeCoords, traction, thickness)
% 2-node boundary element for distributed mechanical traction
%
% nodeCoords : 2 x 2 array [x y]
% traction   : [tx, ty]
% thickness  : thickness
%
% Output:
%   edgeVector : 4 x 1 vector for [ux1 uy1 ux2 uy2]'

x1 = nodeCoords(1,1); y1 = nodeCoords(1,2);
x2 = nodeCoords(2,1); y2 = nodeCoords(2,2);

edgeLength = sqrt((x2-x1)^2 + (y2-y1)^2);

tractionX = traction(1);
tractionY = traction(2);

gaussPoints = [-1/sqrt(3), 1/sqrt(3)];
gaussWeights = [1, 1];

edgeVector = zeros(4,1);

for iGauss = 1:2
    s = gaussPoints(iGauss);

    N = [(1-s)/2, (1+s)/2];
    jacobian1D = edgeLength/2;

    Nmat = [N(1) 0    N(2) 0;
            0    N(1) 0    N(2)];

    edgeVector = edgeVector + Nmat' * [tractionX; tractionY] * jacobian1D * gaussWeights(iGauss) * thickness;
end

end