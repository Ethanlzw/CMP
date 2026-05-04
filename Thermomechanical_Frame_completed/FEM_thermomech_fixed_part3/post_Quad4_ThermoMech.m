function out = post_Quad4_ThermoMech(eleNodesInfo, material, problem, eleSol)
% Postprocess Quad4 thermomechanical element at centroid
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

u = eleSol([1 2 4 5 7 8 10 11]);
T = eleSol([3 6 9 12]);

x = eleNodesInfo(:,2);
y = eleNodesInfo(:,3);

D = getElasticMatrix2D(E, nu, mechType);
Kcond = [kxx 0;
         0   kyy];

s = 0.0;
t = 0.0;

N = 0.25 * [(1+s)*(1+t);
            (1-s)*(1+t);
            (1-s)*(1-t);
            (1+s)*(1-t)];
                
dNds = 0.25 * [(1+t); -(1+t); -(1-t);  (1-t)];
dNdt = 0.25 * [(1+s);  (1-s); -(1-s); -(1+s)];

J = [dNds'*x, dNds'*y;
     dNdt'*x, dNdt'*y];

detJ = det(J);
if detJ <= 0
    error('post_Quad4_ThermoMech found non-positive detJ. Check node ordering or element distortion.');
end

gradN = J \ [dNds'; dNdt'];
dNdx = gradN(1,:)';
dNdy = gradN(2,:)';

Bu = [dNdx(1) 0       dNdx(2) 0       dNdx(3) 0       dNdx(4) 0;
      0       dNdy(1) 0       dNdy(2) 0       dNdy(3) 0       dNdy(4);
      dNdy(1) dNdx(1) dNdy(2) dNdx(2) dNdy(3) dNdx(3) dNdy(4) dNdx(4)];

BT = [dNdx';
      dNdy'];

strain = Bu * u;
Tc = N' * T;
strainTh = alpha * (Tc - Tref) * [1; 1; 0];
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
