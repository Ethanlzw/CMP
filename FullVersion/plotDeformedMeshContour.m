function plotDeformedMeshContour(model, solution, scaleFactor, figTitle)
% Plot displacement magnitude contour on deformed mesh

if nargin < 3 || isempty(scaleFactor)
    scaleFactor = 1.0;
end
if nargin < 4 || isempty(figTitle)
    figTitle = 'Deformed Mesh Contour';
end

if ~(isequal(model.fieldNames, {'ux','uy'}) || isequal(model.fieldNames, {'ux','uy','T'}))
    error('plotDeformedMeshContour requires displacement fields ux and uy.');
end

nodal = extractNodalFields(model, solution);
ux = nodal.ux;
uy = nodal.uy;
displacementMagnitude = sqrt(ux.^2 + uy.^2);

nodeCoords = model.nodesInfo(:,2:3);
deformedCoords = [nodeCoords(:,1) + scaleFactor * ux, ...
                  nodeCoords(:,2) + scaleFactor * uy];

figure('Name', figTitle);
hold on;
axis equal;
box on;

elementTypeTable = model.elesType;
elementSetNames = model.elesSetName;
nElementSets = size(elementTypeTable,1);

for iSet = 1:nElementSets
    nElementNodes = elementTypeTable{iSet,4};
    elementSet = model.elesInfo.(elementSetNames{iSet});

    if isempty(elementSet)
        continue;
    end

    faces = elementSet(:,2:1+nElementNodes);

    patch('Faces', faces, ...
          'Vertices', deformedCoords, ...
          'FaceVertexCData', displacementMagnitude, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'k', ...
          'LineWidth', 0.5);
end

colorbar;
title(figTitle);
xlabel('x');
ylabel('y');

hold off;

end