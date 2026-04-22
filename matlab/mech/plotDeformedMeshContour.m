function plotDeformedMeshContour(model, result, fieldName, scale)
% plotDeformedMeshContour  Deformed mesh with nodal field colour contour.
%
%   fieldName – string matching a field in result.fields, e.g. 'uy'
%   scale     – displacement magnification factor (default 1)

if nargin < 4, scale = 1.0; end

if ~isfield(result.fields, fieldName)
    error('plotDeformedMeshContour: field ''%s'' not found in result.fields.', fieldName);
end

coords     = model.nodesInfo(:,2:3) + scale * [result.fields.ux, result.fields.uy];
fieldValue = result.fields.(fieldName);

figure;  hold on;  axis equal;  box on;

for s = 1:size(model.elesSets, 1)
    setField = model.elesSets{s,1};
    elesSet  = model.elesInfo.(setField);
    for e = 1:size(elesSet, 1)
        en = elesSet(e,:);
        patch(coords(en,1), coords(en,2), fieldValue(en), ...
              'EdgeColor', 'k', 'FaceColor', 'interp');
    end
end

colorbar;
title(sprintf('Deformed Contour: %s  (scale ×%.0e)', fieldName, scale));
xlabel('x [m]');  ylabel('y [m]');

end
