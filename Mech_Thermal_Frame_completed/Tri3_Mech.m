function [eleMatrix, eleVector, info] = Tri3_Mech(eleNodesInfo, material, problem, loads, state)
% Tri3 linear mechanical element
% DoF order per node: [ux, uy]
% 自由度顺序:
%   [ux1 uy1 ux2 uy2 ux3 uy3]'
%
% 输入说明:
%   eleNodesInfo : 3 x 3, 每行为 [nodeID x y]
%   material     : [E, nu]
%   problem      : 结构体, 可含 thickness, mechType
%   loads        : 结构体, 可含 bodyForce = [bx; by]
%   state        : 预留接口, 当前纯力学线性问题中未使用
%
% 输出:
%   eleMatrix : 6 x 6 单元刚度矩阵
%   eleVector : 6 x 1 单元等效载荷向量
%   info      : 额外信息, 供调试和后处理使用

% ----- Geometry -----
x1 = eleNodesInfo(1,2); y1 = eleNodesInfo(1,3);
x2 = eleNodesInfo(2,2); y2 = eleNodesInfo(2,3);
x3 = eleNodesInfo(3,2); y3 = eleNodesInfo(3,3);

% ----- Material / problem -----
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

% 平面应力 / 平面应变本构矩阵
D = getElasticMatrix2D(E, nu, mechType);

%% CORE LOGIC: Tri3 面积、b_i c_i、B 矩阵、K_e、体力一致载荷
% 对 Tri3 单元，位移插值是一阶多项式，所以应变在单元内为常数.
A = [1 x1 y1;
     1 x2 y2;
     1 x3 y3];
detA = det(A);
Area = detA / 2;

% b_i, c_i 是形函数导数常数项.
b1 = y2 - y3;  c1 = x3 - x2;
b2 = y3 - y1;  c2 = x1 - x3;
b3 = y1 - y2;  c3 = x2 - x1;

% B 矩阵将节点位移映射为常应变:
%   strain = B * u_e
B = [b1 0  b2 0  b3 0;
     0  c1 0  c2 0  c3;
     c1 b1 c2 b2 c3 b3] / detA;

% ----- Stiffness matrix -----
% Tri3 为常应变单元, B 与 D 在单元内均为常数:
%   K_e = t * A * B^T * D * B
eleMatrix = thickness * Area * (B' * D * B);

% ----- Consistent body-force vector -----
% 常体力 [bx; by] 的一致载荷:
%   f_e = t * A / 3 * [bx by bx by bx by]^T
eleVector = zeros(6,1);

bodyForce = loads.bodyForce(:);
eleVector = thickness * Area / 3 * [bodyForce(1);
                                    bodyForce(2);
                                    bodyForce(1);
                                    bodyForce(2);
                                    bodyForce(1);
                                    bodyForce(2)];


%% ----- Extra info -----
info.Area = Area;
info.B = B;
info.D = D;
info.type = 'Tri3_Mech';

end