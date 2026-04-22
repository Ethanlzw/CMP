function plotDeformedMesh(model, result, scale)
% plotDeformedMesh  Draw the deformed mesh with displacement magnification.
%
%   scale – displacement magnification factor (default 1)

if nargin < 3, scale = 1.0; end

coords = model.nodesInfo(:,2:3) + scale * [result.fields.ux, result.fields.uy];
figure;
plotElementPatches(model, coords, 'w', 'b');
title(sprintf('Deformed Mesh  (scale ×%.0e)', scale));
xlabel('x [m]');  ylabel('y [m]');

end
