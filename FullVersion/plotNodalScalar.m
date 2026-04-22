function plotNodalScalar(model, nodalValues, figTitle)
% Plot nodal scalar field on original mesh

if nargin < 3 || isempty(figTitle)
    figTitle = 'Nodal Scalar';
end

nNodes = size(model.nodesInfo,1);
nodalValues = nodalValues(:);

if numel(nodalValues) ~= nNodes
    error('nodalValues size must match number of nodes.');
end

nodeCoords = model.nodesInfo(:,2:3);

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
          'Vertices', nodeCoords, ...
          'FaceVertexCData', nodalValues, ...
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