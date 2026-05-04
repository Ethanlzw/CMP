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

L = sqrt((x2-x1)^2 + (y2-y1)^2);

tx = traction(1);
ty = traction(2);

gp = [-1/sqrt(3), 1/sqrt(3)];
w  = [1, 1];

edgeVector = zeros(4,1);

for igp = 1:2
    s = gp(igp);

    N = [(1-s)/2, (1+s)/2];
    J = L/2;

    Nmat = [N(1) 0    N(2) 0;
            0    N(1) 0    N(2)];

    edgeVector = edgeVector + Nmat' * [tx; ty] * J * w(igp) * thickness;
end

end