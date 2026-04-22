function [bcMatrix, bcVector] = assembleNaturalBC_Thermal(model)
% Assemble thermal natural BC contributions
%
% model.natBC rows:
%   { [node1,node2,...], edgeDataStruct }
%
% edgeDataStruct.type = 'flux' or 'convection'

nNodeDofs = model.elesType{1,3};
if nNodeDofs ~= 1
    error('assembleNaturalBC_Thermal is only for thermal models with 1 DoF per node.');
end

nNodes = size(model.nodesInfo,1);
bcMatrix = zeros(nNodes,nNodes);
bcVector = zeros(nNodes,1);

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

        [edgeMatrix, edgeVector] = Edge2_ThermalFlux(nodeCoords, edgeData, thickness);

        edgeDofs = [nodeA; nodeB];
        bcMatrix(edgeDofs,edgeDofs) = bcMatrix(edgeDofs,edgeDofs) + edgeMatrix;
        bcVector(edgeDofs) = bcVector(edgeDofs) + edgeVector;
    end
end

end