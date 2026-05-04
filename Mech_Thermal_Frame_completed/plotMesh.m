function plotMesh(model, figTitleStr, showNodeIDs, showElementIDs)
% Plot mixed Tri3/Quad4 mesh

if nargin < 2
    figTitleStr = 'Mesh';
end
if nargin < 3
    showNodeIDs = false;
end
if nargin < 4
    showElementIDs = false;
end

figure;
hold on;
axis equal;
box on;
view(2);

elesType = model.elesType;
elesSetName = model.elesSetName;
elesSetNum = size(elesType,1);

for setI = 1:elesSetNum
    eleName = elesType{setI,2};
    nEleNodes = elesType{setI,4};
    elesSet = model.elesInfo.(elesSetName{setI});
    nSetEles = size(elesSet,1);

    for e = 1:nSetEles
        eleID = elesSet(e,1);
        eleNodeIDs = elesSet(e,2:1+nEleNodes);

        coords = zeros(nEleNodes,2);
        for a = 1:nEleNodes
            idx = find(model.nodesInfo(:,1)==eleNodeIDs(a),1);
            coords(a,:) = model.nodesInfo(idx,2:3);
        end

        patch(coords(:,1), coords(:,2), 'w', ...
              'EdgeColor', 'k', ...
              'LineWidth', 1.0, ...
              'FaceColor', 'none');

        if showElementIDs
            c = mean(coords,1);
            text(c(1), c(2), sprintf('%d', eleID), ...
                'Color', 'r', 'HorizontalAlignment', 'center');
        end
    end
end

if showNodeIDs
    for n = 1:size(model.nodesInfo,1)
        text(model.nodesInfo(n,2), model.nodesInfo(n,3), ...
            sprintf('%d', model.nodesInfo(n,1)), ...
            'Color', 'b', ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'right');
    end
end

title(figTitleStr);
xlabel('x');
ylabel('y');

end