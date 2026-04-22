function plotMesh(model, figTitle, showNodeID, showEleID)
% Plot FE mesh

if nargin < 2 || isempty(figTitle)
    figTitle = 'FE Mesh';
end
if nargin < 3 || isempty(showNodeID)
    showNodeID = false;
end
if nargin < 4 || isempty(showEleID)
    showEleID = false;
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
          'FaceColor', 'none', ...
          'EdgeColor', 'k', ...
          'LineWidth', 1.0);
end

if showNodeID
    for iNode = 1:size(model.nodesInfo,1)
        text(nodeCoords(iNode,1), nodeCoords(iNode,2), sprintf('%d', model.nodesInfo(iNode,1)), ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left');
    end
end

if showEleID
    for iSet = 1:nElementSets
        nElementNodes = elementTypeTable{iSet,4};
        elementSet = model.elesInfo.(elementSetNames{iSet});

        for iElement = 1:size(elementSet,1)
            elementID = elementSet(iElement,1);
            elementNodeIDs = elementSet(iElement,2:1+nElementNodes);
            xy = model.nodesInfo(elementNodeIDs,2:3);
            center = mean(xy,1);

            text(center(1), center(2), sprintf('%d', elementID), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle');
        end
    end
end

title(figTitle);
xlabel('x');
ylabel('y');

hold off;

end