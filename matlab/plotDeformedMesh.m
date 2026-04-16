function plotDeformedMesh(model, result, scale)
% Plot deformed mesh

if nargin < 3
    scale = 1.0;
end

ux = result.fields.ux;
uy = result.fields.uy;

nodesDef = model.nodesInfo(:,2:3) + scale * [ux, uy];

figure;
hold on;
axis equal;
box on;

for setI = 1:numel(model.elesSetName)
    setName = model.elesSetName{setI};
    elesSet = model.elesInfo.(setName);

    for eleI = 1:size(elesSet,1)
        eleNodes = elesSet(eleI,:);
        xy = nodesDef(eleNodes, :);
        patch(xy(:,1), xy(:,2), 'w', 'EdgeColor', 'b');
    end
end

title('Deformed Mesh');
xlabel('x');
ylabel('y');

end