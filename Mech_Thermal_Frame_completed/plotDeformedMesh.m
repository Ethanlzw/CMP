function plotDeformedMesh(model, solution, scaleFactor, figTitleStr)
% Plot undeformed and deformed mixed mesh
%
% scaleFactor: deformation magnification factor

if nargin < 3 || isempty(scaleFactor)
    scaleFactor = 1.0;
end

if nargin < 4
    figTitleStr = 'Deformed Mesh';
end

nNodeDoF = model.elesType{1,3};
if nNodeDoF < 2
    error('Deformed mesh plot requires mechanical DoFs.');
end

numNodes = size(model.nodesInfo,1);
coords0 = model.nodesInfo(:,2:3);
coordsDef = coords0;

for n = 1:numNodes
    base = (n-1)*nNodeDoF;
    ux = solution(base+1);
    uy = solution(base+2);
    coordsDef(n,1) = coordsDef(n,1) + scaleFactor*ux;
    coordsDef(n,2) = coordsDef(n,2) + scaleFactor*uy;
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
        eleNodeIDs = elesSet(e,2:1+nEleNodes);

        face0 = coords0(eleNodeIDs, :);
        faceD = coordsDef(eleNodeIDs, :);

        % undeformed
        patch(face0(:,1), face0(:,2), 'w', ...
              'EdgeColor', [0.6 0.6 0.6], ...
              'LineStyle', '--', ...
              'LineWidth', 0.8, ...
              'FaceColor', 'none');

        % deformed
        patch(faceD(:,1), faceD(:,2), 'w', ...
              'EdgeColor', 'k', ...
              'LineWidth', 1.2, ...
              'FaceColor', 'none');
    end
end

title(figTitleStr);
xlabel('x');
ylabel('y');

legend({'Undeformed','Deformed'});
end