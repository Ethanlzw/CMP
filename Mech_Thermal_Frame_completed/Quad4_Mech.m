  function [eleMatrix, eleVector, info] = Quad4_Mech(eleNodesInfo, material, problem, loads, state)
% Quad4 linear mechanical element
% DoF order per node: [ux, uy]
% 自由度顺序:
%   [ux1 uy1 ux2 uy2 ux3 uy3 ux4 uy4]'
% 说明:
%   采用 2x2 Gauss 积分计算刚度和一致体力载荷.
%   结点顺序假定为标准逆时针:
%       1(-1,-1), 2(1,-1), 3(1,1), 4(-1,1)

x = eleNodesInfo(:,2);
y = eleNodesInfo(:,3);

E  = material(1);
nu = material(2);

if isfield(problem,'thickness')
    thickness = problem.thickness;
else
    thickness = 1.0;
end

if isfield(problem,'mechType')
    mechType = problem.mechType;
else
    mechType = 'planeStress';
end

bodyForce = [0; 0];
if nargin >= 4 && ~isempty(loads) && isfield(loads,'bodyForce')
    bodyForce = loads.bodyForce(:);
    if numel(bodyForce) ~= 2
        error('Mechanical bodyForce must be a 2-component vector [bx; by].');
    end
end

D = getElasticMatrix2D(E, nu, mechType);  

%% CORE LOGIC: Q4 的 N, dN/ds, dN/dt, Jacobian, B 矩阵, 2x2 Gauss 积分
% 2x2 Gauss points
gaussPts = [-1/sqrt(3), -1/sqrt(3);
            -1/sqrt(3),  1/sqrt(3);
             1/sqrt(3), -1/sqrt(3);
             1/sqrt(3),  1/sqrt(3)];
weights = [1; 1; 1; 1];

eleMatrix = zeros(8,8);
eleVector = zeros(8,1);

%% ----- Calculate eleMatrix and eleVector using for loop -----

for igp = 1:4
    s = gaussPts(igp, 1);
    t = gaussPts(igp, 2);
    w = weights(igp);
    
    % 形函数
    N = 0.25 * [(1+s)*(1+t);
                (1-s)*(1+t);
                (1-s)*(1-t);
                (1+s)*(1-t)];

    dNds = 0.25 * [(1+t); -(1+t); -(1-t);  (1-t)];
    dNdt = 0.25 * [(1+s);  (1-s); -(1-s); -(1+s)];

    % Jacobian
    J = [dNds'*x, dNds'*y;
         dNdt'*x, dNdt'*y];
    detJ = det(J);

    % 将形函数导数从 (s,t) 映射到 (x,y)
    gradN = J \ [dNds; dNdt];
    dNdx = gradN(1, :);
    dNdy = gradN(2, :);

    % B 矩阵: strain = B * u_e
    B = [dNdx(1)      0    dNdx(2)   0      dNdx(3)     0      dNdx(4)    0;
           0      dNdy(1)    0     dNdy(2)    0      dNdy(3)    0       dNdy(4);
         dNdy(1)  dNdx(1)  dNdy(2) dNdx(2)  dNdy(3)  dNdx(3)   dNdy(4)  dNdx(4)];

    % 单元刚度积分
    eleMatrix = eleMatrix + (B' * D * B) * detJ * w * thickness;

    % 体力一致载荷积分
    Nmat = [N(1)   0    N(2)    0    N(3)   0    N(4)   0;
             0    N(1)    0    N(2)   0    N(3)   0    N(4)];

    eleVector = eleVector + Nmat' * bodyForce * detJ * w * thickness;   
end

info.D = D;
info.type = 'Quad4_Mech';

end
