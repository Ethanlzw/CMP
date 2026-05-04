function plotNodalScalar(model, nodalValues, figTitleStr)
% Plot nodal scalar field on mixed Tri3/Quad4 mesh

figure;
hold on;
axis equal;
box on;
view(2);

nodesInfo = model.nodesInfo;
elesType = model.elesType;
elesSetName = model.elesSetName;

x = nodesInfo(:,2);
y = nodesInfo(:,3);

elesSetNum = size(elesType,1);

for setI = 1:elesSetNum
    eleName = elesType{setI,2};
    nEleNodes = elesType{setI,4};
    elesSet = model.elesInfo.(elesSetName{setI});

    nSetEles = size(elesSet,1);
    faces = zeros(nSetEles, nEleNodes);

    for e = 1:nSetEles
        eleNodeIDs = elesSet(e,2:1+nEleNodes);
        faces(e,:) = getElementFaces2D(eleNodeIDs, eleName);
    end

    patch('Faces', faces, ...
          'Vertices', [x y], ...
          'FaceVertexCData', nodalValues, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'k', ...
          'LineWidth', 0.8);
end

colorbar;
title(figTitleStr);
xlabel('x');
ylabel('y');

end