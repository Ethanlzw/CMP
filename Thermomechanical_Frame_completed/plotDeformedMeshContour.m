function plotDeformedMeshContour(model, solution, scaleFactor, figTitleStr)
% Plot deformed mixed mesh with nodal displacement magnitude contour

if nargin < 3 || isempty(scaleFactor)
    scaleFactor = 1.0;
end

if nargin < 4
    figTitleStr = 'Deformed Mesh with Displacement Magnitude';
end

nNodeDoF = model.elesType{1,3};
if nNodeDoF < 2
    error('Deformed mesh contour requires mechanical DoFs.');
end

numNodes = size(model.nodesInfo,1);
coordsDef = model.nodesInfo(:,2:3);
uMag = zeros(numNodes,1);

for n = 1:numNodes
    base = (n-1)*nNodeDoF;
    ux = solution(base+1);
    uy = solution(base+2);
    coordsDef(n,1) = coordsDef(n,1) + scaleFactor*ux;
    coordsDef(n,2) = coordsDef(n,2) + scaleFactor*uy;
    uMag(n) = sqrt(ux^2 + uy^2);
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
    faces = zeros(nSetEles, nEleNodes);

    for e = 1:nSetEles
        eleNodeIDs = elesSet(e,2:1+nEleNodes);
        faces(e,:) = getElementFaces2D(eleNodeIDs, eleName);
    end

    patch('Faces', faces, ...
          'Vertices', coordsDef, ...
          'FaceVertexCData', uMag, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'k', ...
          'LineWidth', 0.8);
end

colorbar;
title(figTitleStr);
xlabel('x');
ylabel('y');

end