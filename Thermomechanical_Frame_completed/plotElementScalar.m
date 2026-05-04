function plotElementScalar(model, elementResults, fieldName, figTitleStr)
% Plot element-wise scalar field on mixed Tri3/Quad4 mesh

figure;
hold on;
axis equal;
box on;
view(2);

vals = zeros(numel(elementResults),1);
for e = 1:numel(elementResults)
    vals(e) = elementResults(e).value.(fieldName);
end

for e = 1:numel(elementResults)
    eleNodeIDs = elementResults(e).nodeIDs;
    val = elementResults(e).value.(fieldName);

    coords = zeros(numel(eleNodeIDs),2);
    for a = 1:numel(eleNodeIDs)
        idx = find(model.nodesInfo(:,1)==eleNodeIDs(a),1);
        coords(a,:) = model.nodesInfo(idx,2:3);
    end

    patch(coords(:,1), coords(:,2), val, ...
          'EdgeColor', 'k', ...
          'LineWidth', 0.8);
end

colorbar;
caxis([min(vals), max(vals)]);
title(figTitleStr);
xlabel('x');
ylabel('y');

end