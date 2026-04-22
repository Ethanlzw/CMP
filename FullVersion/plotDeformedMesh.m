function plotDeformedMesh(model, solution, scaleFactor, figTitle)
% Plot undeformed and deformed mesh for mechanical or thermo-mechanical models

if nargin < 3 || isempty(scaleFactor)
    scaleFactor = 1.0;
end
if nargin < 4 || isempty(figTitle)
    figTitle = 'Deformed Mesh';
end

if ~(isequal(model.fieldNames, {'ux','uy'}) || isequal(model.fieldNames, {'ux','uy','T'}))
    error('plotDeformedMesh requires displacement fields ux and uy.');
end

nodal = extractNodalFields(model, solution);
ux = nodal.ux;
uy = nodal.uy;

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
          'Vertices', nodeCoords, ...
          'FaceColor', 'none', ...
          'EdgeColor', [0.6 0.6 0.6], ...
          'LineWidth', 0.75);

    patch('Faces', faces, ...
          'Vertices', deformedCoords, ...
          'FaceColor', 'none', ...
          'EdgeColor', 'b', ...
          'LineWidth', 1.0);
end

title(sprintf('%s (scale = %.3g)', figTitle, scaleFactor));
xlabel('x');
ylabel('y');

hold off;

end