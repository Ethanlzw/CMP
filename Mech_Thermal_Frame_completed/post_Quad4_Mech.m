function out = post_Quad4_Mech(eleNodesInfo, material, problem, eleSol)
% Postprocess Quad4 mechanical element at centroid
% eleSol = [ux1 uy1 ux2 uy2 ux3 uy3 ux4 uy4]'
% 说明:
%   Quad4 单元的应变在单元内随位置变化.
%   这里取单元中心 (s,t) = (0,0) 处的结果作为代表值.
%   这种做法简单直观，适合教学和基础后处理.

x = eleNodesInfo(:,2);
y = eleNodesInfo(:,3);

E  = material(1);
nu = material(2);

if isfield(problem,'mechType')
    mechType = problem.mechType;
else
    mechType = 'planeStress';
end

D = getElasticMatrix2D(E, nu, mechType);

%% MODIFICATION: 由单元位移求 strain, 再由 stress = D * strain

% 单元中心
s = 0.0; 
t = 0.0;

% 形函数对自然坐标的导数
dNds = 0.25 * [(1+t); -(1+t); -(1-t);  (1-t)];
dNdt = 0.25 * [(1+s);  (1-s); -(1-s); -(1+s)];

% Jacobian
J = [dNds'*x, dNds'*y;
     dNdt'*x, dNdt'*y];

% 从 (s,t) 导数映射到 (x,y) 导数
gradN = inv(J) * [dNds'; dNdt'];
dNdx = gradN(1,:)';
dNdy = gradN(2,:)';

B = [dNdx(1)   0     dNdx(2)   0     dNdx(3)   0     dNdx(4)   0;
       0     dNdy(1)   0     dNdy(2)   0     dNdy(3)   0     dNdy(4);
     dNdy(1) dNdx(1) dNdy(2) dNdx(2) dNdy(3) dNdx(3) dNdy(4) dNdx(4)];

strain = B * eleSol;
stress = D * strain;

out.centroid = mean(eleNodesInfo(:,2:3),1);
out.strain = strain;
out.stress = stress;
out.vonMises = sqrt(stress(1)^2 - stress(1)*stress(2) + stress(2)^2 + 3*stress(3)^2);

end
