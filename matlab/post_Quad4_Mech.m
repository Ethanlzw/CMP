function out = post_Quad4_Mech(eleNodesInfo, eleU, material, problem, state)
% Post-process Quad4 mechanical element
% Usually evaluate at element center

x = eleNodesInfo(:,2);
y = eleNodesInfo(:,3);

E  = material(1);
nu = material(2);

if isfield(problem, 'mechType')
    mechType = problem.mechType;
else
    mechType = 'planeStress';
end

D = getElasticMatrix2D(E, nu, mechType);

s = 0.0;
t = 0.0;

dNds = 0.25 * [-(1-t),  (1-t), (1+t), -(1+t)];
dNdt = 0.25 * [-(1-s), -(1+s), (1+s),  (1-s)];

J = [dNds; dNdt] * [x y];
dNdxy = J \ [dNds; dNdt];

B = zeros(3,8);
for i = 1:4
    B(:,2*i-1:2*i) = [dNdxy(1,i), 0;
                      0, dNdxy(2,i);
                      dNdxy(2,i), dNdxy(1,i)];
end

strain = B * eleU;
stress = D * strain;

out = struct();
out.strain = strain;
out.stress = stress;

end