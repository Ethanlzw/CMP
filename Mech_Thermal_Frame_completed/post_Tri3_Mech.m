function out = post_Tri3_Mech(eleNodesInfo, material, problem, eleSol)
% Postprocess Tri3 mechanical element
% eleSol = [ux1 uy1 ux2 uy2 ux3 uy3]'

x1 = eleNodesInfo(1,2); y1 = eleNodesInfo(1,3);
x2 = eleNodesInfo(2,2); y2 = eleNodesInfo(2,3);
x3 = eleNodesInfo(3,2); y3 = eleNodesInfo(3,3);

E  = material(1);
nu = material(2);

if isfield(problem,'mechType')
    mechType = problem.mechType;
else
    mechType = 'planeStress';
end

D = getElasticMatrix2D(E, nu, mechType);

A = [1 x1 y1;
     1 x2 y2;
     1 x3 y3];
detA = det(A);

b1 = y2-y3; b2 = y3-y1; b3 = y1-y2;
c1 = x3-x2; c2 = x1-x3; c3 = x2-x1;

B = [b1 0 b2 0 b3 0;
     0 c1 0 c2 0 c3;
     c1 b1 c2 b2 c3 b3] / detA;

strain = B * eleSol;
stress = D * strain;

out.centroid = mean(eleNodesInfo(:,2:3),1);
out.strain = strain;
out.stress = stress;
out.vonMises = sqrt(stress(1)^2 - stress(1)*stress(2) + stress(2)^2 + 3*stress(3)^2);

end