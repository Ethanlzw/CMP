function plotDeformedMeshContour(model, result, fieldName, scale)
% Plot deformed mesh with nodal contour

if nargin < 4
    scale = 1.0;
end

ux = result.fields.ux;
uy = result.fields.uy;
fieldValue = result.fields.(fieldName);

nodesDef = model.nodesInfo(:,2:3) + scale * [ux, uy];

figure;
hold on;
axis equal;
box on;

for setI = 1:numel(model.elesSetName)
    setName = model.elesSetName{setI};
    elesSet = model.elesInfo.(setName);

    for eleI = 1:size(elesSet,1)
        eleNodes = elesSet(eleI,:);
        xy = nodesDef(eleNodes, :);
        c  = fieldValue(eleNodes);

        patch(xy(:,1), xy(:,2), c, 'EdgeColor', 'k', 'FaceColor', 'interp');
    end
end

colorbar;
title(['Deformed Mesh Contour: ' fieldName]);
xlabel('x');
ylabel('y');

end