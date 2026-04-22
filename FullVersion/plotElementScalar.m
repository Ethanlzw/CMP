function plotElementScalar(model, elementValues, figTitle)
% Plot element-wise scalar field on original mesh

if nargin < 3 || isempty(figTitle)
    figTitle = 'Element Scalar';
end

elementValues = elementValues(:);

elementTypeTable = model.elesType;
elementSetNames = model.elesSetName;
nElementSets = size(elementTypeTable,1);

nElementsTotal = 0;
for iSet = 1:nElementSets
    elementSet = model.elesInfo.(elementSetNames{iSet});
    nElementsTotal = nElementsTotal + size(elementSet,1);
end

if numel(elementValues) ~= nElementsTotal
    error('elementValues size must match number of elements.');
end

nodeCoords = model.nodesInfo(:,2:3);

figure('Name', figTitle);
hold on;
axis equal;
box on;

iValue = 0;
for iSet = 1:nElementSets
    nElementNodes = elementTypeTable{iSet,4};
    elementSet = model.elesInfo.(elementSetNames{iSet});
    nElementsThisSet = size(elementSet,1);

    if isempty(elementSet)
        continue;
    end

    faces = elementSet(:,2:1+nElementNodes);
    valuesThisSet = elementValues(iValue+1 : iValue+nElementsThisSet);

    patch('Faces', faces, ...
          'Vertices', nodeCoords, ...
          'FaceVertexCData', valuesThisSet, ...
          'FaceColor', 'flat', ...
          'EdgeColor', 'k', ...
          'LineWidth', 0.5);

    iValue = iValue + nElementsThisSet;
end

colorbar;
title(figTitle);
xlabel('x');
ylabel('y');

hold off;

end