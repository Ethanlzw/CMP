function out = post_Quad4_Thermal(eleNodesInfo, eleU, material, problem, state)
% Post-process Quad4 thermal element at center

x = eleNodesInfo(:,2);
y = eleNodesInfo(:,3);

k = material(1);

s = 0.0;
t = 0.0;

dNds = 0.25 * [-(1-t),  (1-t), (1+t), -(1+t)];
dNdt = 0.25 * [-(1-s), -(1+s), (1+s),  (1-s)];

J = [dNds; dNdt] * [x y];
dNdxy = J \ [dNds; dNdt];

B = dNdxy;

gradT = B * eleU;
heatFlux = -k * gradT;

out = struct();
out.gradT = gradT;
out.heatFlux = heatFlux;

end