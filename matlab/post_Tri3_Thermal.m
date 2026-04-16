function out = post_Tri3_Thermal(eleNodesInfo, eleU, material, problem, state)
% Post-process Tri3 thermal element

x1 = eleNodesInfo(1,2); y1 = eleNodesInfo(1,3);
x2 = eleNodesInfo(2,2); y2 = eleNodesInfo(2,3);
x3 = eleNodesInfo(3,2); y3 = eleNodesInfo(3,3);

k = material(1);

A2 = det([1 x1 y1;
          1 x2 y2;
          1 x3 y3]);

A = A2 / 2;

b1 = y2 - y3;  c1 = x3 - x2;
b2 = y3 - y1;  c2 = x1 - x3;
b3 = y1 - y2;  c3 = x2 - x1;

B = 1 / (2*A) * [b1 b2 b3;
                 c1 c2 c3];

gradT = B * eleU;
heatFlux = -k * gradT;

out = struct();
out.gradT = gradT;
out.heatFlux = heatFlux;

end