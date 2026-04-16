function out = post_Tri3_Mech(eleNodesInfo, eleU, material, problem, state)
% Post-process Tri3 mechanical element

[~, ~, info] = Tri3_Mech(eleNodesInfo, material, problem, struct(), state);

B = info.B;
D = info.D;

strain = B * eleU;
stress = D * strain;

out = struct();
out.strain = strain;
out.stress = stress;

end