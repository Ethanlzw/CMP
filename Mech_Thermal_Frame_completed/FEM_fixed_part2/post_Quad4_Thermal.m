function out = post_Quad4_Thermal(eleNodesInfo, material, problem, eleSol)
% Postprocess Quad4 thermal element at centroid
% eleSol = [T1 T2 T3 T4]'
% Local node order:
%   1(-1,-1), 2(1,-1), 3(1,1), 4(-1,1)

x = eleNodesInfo(:,2);
y = eleNodesInfo(:,3);

if numel(material) == 1
    kxx = material(1);
    kyy = material(1);
else
    kxx = material(1);
    kyy = material(2);
end

Kcond = [kxx 0;
         0   kyy];

s = 0.0;
t = 0.0;

N = 0.25 * [(1-s)*(1-t);
            (1+s)*(1-t);
            (1+s)*(1+t);
            (1-s)*(1+t)];

dNds = 0.25 * [-(1-t),  (1-t),  (1+t), -(1+t)];
dNdt = 0.25 * [-(1-s), -(1+s),  (1+s),  (1-s)];

J = [dNds*x, dNds*y;
     dNdt*x, dNdt*y];
detJ = det(J);

if detJ <= 0
    error('post_Quad4_Thermal detected non-positive detJ = %.6e. Check node order or element distortion.', detJ);
end

gradN = J \ [dNds; dNdt];
gradT = gradN * eleSol;
q = -Kcond * gradT;

out.centroid = mean(eleNodesInfo(:,2:3),1);
out.T = N' * eleSol;
out.gradT = gradT;
out.heatFlux = q;

end
