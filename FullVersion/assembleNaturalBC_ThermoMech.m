function [bcMatrix, bcVector] = assembleNaturalBC_ThermoMech(model)
% Assemble thermo-mechanical natural BC contributions
%
% model.natBC rows:
%   { [node1,node2,...], edgeDataStruct }
%
% edgeDataStruct may include:
%   .traction    = [tx; ty]
%   .thermalType = 'flux' or 'convection'
%   .qn / .h / .Tinf

nNodeDoFs = model.elesType{1,3};
if nNodeDoFs ~= 3
    error('assembleNaturalBC_ThermoMech is only for ux-uy-T models.');
end

nNodes = size(model.nodesInfo,1);
nTotalDoFs = nNodes * nNodeDoFs;

bcMatrix = zeros(nTotalDoFs,nTotalDoFs);
bcVector = zeros(nTotalDoFs,1);

if ~isfield(model,'natBC') || isempty(model.natBC)
    return;
end

if isfield(model,'problem') && isfield(model.problem,'thickness')
    thickness = model.problem.thickness;
else
    thickness = 1.0;
end

for iBC = 1:size(model.natBC,1)
    bcNodeIDs = model.natBC{iBC,1};
    edgeData = model.natBC{iBC,2};

    for iEdge = 1:(numel(bcNodeIDs)-1)
        nodeA = bcNodeIDs(iEdge);
        nodeB = bcNodeIDs(iEdge+1);

        edgeInfo = findEdgeOwner(model, nodeA, nodeB);
        if ~edgeInfo.found
            error('Boundary edge (%d,%d) was not found in any element.', nodeA, nodeB);
        end

        idxA = find(model.nodesInfo(:,1) == nodeA, 1);
        idxB = find(model.nodesInfo(:,1) == nodeB, 1);

        nodeCoords = [model.nodesInfo(idxA,2:3);
                      model.nodesInfo(idxB,2:3)];

        [edgeMatrix, edgeVector] = Edge2_ThermoMech(nodeCoords, edgeData, thickness);

        edgeDoFs = [(nodeA-1)*nNodeDoFs + 1;
                    (nodeA-1)*nNodeDoFs + 2;
                    (nodeA-1)*nNodeDoFs + 3;
                    (nodeB-1)*nNodeDoFs + 1;
                    (nodeB-1)*nNodeDoFs + 2;
                    (nodeB-1)*nNodeDoFs + 3];

        bcMatrix(edgeDoFs,edgeDoFs) = bcMatrix(edgeDoFs,edgeDoFs) + edgeMatrix;
        bcVector(edgeDoFs) = bcVector(edgeDoFs) + edgeVector;
    end
end

end