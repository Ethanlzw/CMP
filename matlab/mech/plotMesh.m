function plotMesh(model)
% plotMesh  Draw the undeformed mesh.

coords = model.nodesInfo(:, 2:3);
figure;
plotElementPatches(model, coords, 'w', 'k');
title('Undeformed Mesh');
xlabel('x [m]');  ylabel('y [m]');

end
