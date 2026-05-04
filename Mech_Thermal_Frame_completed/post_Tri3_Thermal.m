function out = post_Tri3_Thermal(eleNodesInfo, material, problem, eleSol)
% Postprocess Tri3 thermal element
% eleSol = [T1 T2 T3]'

x1 = eleNodesInfo(1,2); y1 = eleNodesInfo(1,3);
x2 = eleNodesInfo(2,2); y2 = eleNodesInfo(2,3);
x3 = eleNodesInfo(3,2); y3 = eleNodesInfo(3,3);

kxx = material(1);
kyy = material(2);

A = [1 x1 y1;
     1 x2 y2;
     1 x3 y3];
detA = det(A);

b1 = y2-y3; b2 = y3-y1; b3 = y1-y2;
c1 = x3-x2; c2 = x1-x3; c3 = x2-x1;

B = [b1 b2 b3;
     c1 c2 c3] / detA;

Kcond = [kxx 0;
         0   kyy];

gradT = B * eleSol;
q = -Kcond * gradT;

out.centroid = mean(eleNodesInfo(:,2:3),1);
out.T = mean(eleSol);
out.gradT = gradT;
out.heatFlux = q;

end