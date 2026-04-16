function plotMesh(model)
% Plot undeformed mesh

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
        patch(xy(:,1), xy(:,2), 'w', 'EdgeColor', 'k');
    end
end

title('Mesh');
xlabel('x');
ylabel('y');

end