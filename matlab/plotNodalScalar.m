function plotNodalScalar(model, result, fieldName)
% Plot nodal scalar contour

fieldValue = result.fields.(fieldName);

figure;
hold on;
axis equal;
box on;

for setI = 1:numel(model.elesSetName)
    setName = model.elesSetName{setI};
    elesSet = model.elesInfo.(setName);

    for eleI = 1:size(elesSet,1)
        eleNodes = elesSet(eleI,:);
        xy = model.nodesInfo(eleNodes, 2:3);
        c = fieldValue(eleNodes);

        patch(xy(:,1), xy(:,2), c, 'EdgeColor', 'k', 'FaceColor', 'interp');
    end
end

colorbar;
title(['Nodal Scalar: ' fieldName]);
xlabel('x');
ylabel('y');

end