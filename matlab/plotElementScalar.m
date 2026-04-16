function plotElementScalar(model, eleScalar)
% Plot element-wise scalar

figure;
hold on;
axis equal;
box on;

for setI = 1:numel(model.elesSetName)
    setName = model.elesSetName{setI};
    elesSet = model.elesInfo.(setName);

    if ~isfield(eleScalar, setName)
        continue;
    end

    vals = eleScalar.(setName);

    for eleI = 1:size(elesSet,1)
        eleNodes = elesSet(eleI,:);
        xy = model.nodesInfo(eleNodes, 2:3);
        patch(xy(:,1), xy(:,2), vals(eleI), 'EdgeColor', 'k');
    end
end

colorbar;
title('Element Scalar');

end