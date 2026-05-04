function out = post_Tri3_ThermoMech(eleNodesInfo, material, problem, eleSol)
% Postprocess Tri3 thermomechanical element
% eleSol order per node: [ux uy T]

E     = material(1);
nu    = material(2);
alpha = material(3);
kxx   = material(4);
kyy   = material(5);

if isfield(problem,'mechType')
    mechType = problem.mechType;
else
    mechType = 'planeStress';
end

if isfield(problem,'Tref')
    Tref = problem.Tref;
else
    Tref = 0.0;
end

u = eleSol([1 2 4 5 7 8]);
T = eleSol([3 6 9]);

x1 = eleNodesInfo(1,2); y1 = eleNodesInfo(1,3);
x2 = eleNodesInfo(2,2); y2 = eleNodesInfo(2,3);
x3 = eleNodesInfo(3,2); y3 = eleNodesInfo(3,3);

D = getElasticMatrix2D(E, nu, mechType);
Kcond = [kxx 0;
         0   kyy];

A = [1 x1 y1;
     1 x2 y2;
     1 x3 y3];
detA = det(A);

b1 = y2-y3; b2 = y3-y1; b3 = y1-y2;
c1 = x3-x2; c2 = x1-x3; c3 = x2-x1;

Bu = [b1 0 b2 0 b3 0;
      0 c1 0 c2 0 c3;
      c1 b1 c2 b2 c3 b3] / detA;

BT = [b1 b2 b3;
      c1 c2 c3] / detA;

strain = Bu * u;
Tc = mean(T);
strainTh = alpha * (Tc - Tref) * [1;1;0];
stress = D * (strain - strainTh);

gradT = BT * T;
q = -Kcond * gradT;

out.centroid = mean(eleNodesInfo(:,2:3),1);
out.u = u;
out.T = Tc;
out.strain = strain;
out.thermalStrain = strainTh;
out.stress = stress;
out.vonMises = sqrt(stress(1)^2 - stress(1)*stress(2) + stress(2)^2 + 3*stress(3)^2);
out.gradT = gradT;
out.heatFlux = q;

end