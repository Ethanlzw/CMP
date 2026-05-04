function out = post_Quad4_Thermal(eleNodesInfo, material, problem, eleSol)
% Postprocess Quad4 thermal element at centroid
% eleSol = [T1 T2 T3 T4]'

x = eleNodesInfo(:,2);
y = eleNodesInfo(:,3);

kxx = material(1);
kyy = material(2);

Kcond = [kxx 0;
         0   kyy];

s = 0.0; t = 0.0;

N = 0.25 * [(1+s)*(1+t);
                (1-s)*(1+t);
                (1-s)*(1-t);
                (1+s)*(1-t)];

dNds = 0.25 * [(1+t); -(1+t); -(1-t);  (1-t)];
dNdt = 0.25 * [(1+s);  (1-s); -(1-s); -(1+s)];

J = [dNds'*x, dNds'*y;
     dNdt'*x, dNdt'*y];

gradN = inv(J) * [dNds'; dNdt'];
gradT = gradN * eleSol;
q = -Kcond * gradT;

out.centroid = mean(eleNodesInfo(:,2:3),1);
out.T = N' * eleSol;
out.gradT = gradT;
out.heatFlux = q;

end
