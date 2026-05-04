function [bcMatrix, bcVector] = assembleNaturalBC_Thermal(model)
% Assemble thermal natural BC contributions from model.natBC
%
% model.natBC rows:
%   { [node1,node2,...], [M,S] }

nNodeDoF = model.elesType{1,3};
if nNodeDoF ~= 1
    error('assembleNaturalBC_Thermal is only for thermal models with 1 DoF per node.');
end

numNodes = size(model.nodesInfo,1);
bcMatrix = zeros(numNodes*nNodeDoF);
bcVector = zeros(numNodes*nNodeDoF,1);

if ~isfield(model,'natBC') || isempty(model.natBC)
    return;
end

if isfield(model,'problem') && isfield(model.problem,'thickness')
    thickness = model.problem.thickness;
else
    thickness = 1.0;
end

for ibc = 1:size(model.natBC,1)
    bcNodes = model.natBC{ibc,1};
    bcVals  = model.natBC{ibc,2};

    for i = 1:(numel(bcNodes)-1)
        nodeA = bcNodes(i);
        nodeB = bcNodes(i+1);

        edgeInfo = findEdgeOwner(model, nodeA, nodeB);
        if ~edgeInfo.found
            error('Boundary edge (%d,%d) was not found in any element.', nodeA, nodeB);
        end

        idxA = find(model.nodesInfo(:,1)==nodeA,1);
        idxB = find(model.nodesInfo(:,1)==nodeB,1);
        nodeCoords = [model.nodesInfo(idxA,2:3);
                      model.nodesInfo(idxB,2:3)];

        [edgeMatrix, edgeVector] = Edge2_Thermal(nodeCoords, bcVals, thickness);

        edgeDoFs = [(nodeA-1)*nNodeDoF+1;
                    (nodeB-1)*nNodeDoF+1];

        bcMatrix(edgeDoFs, edgeDoFs) = bcMatrix(edgeDoFs, edgeDoFs) + edgeMatrix;
        bcVector(edgeDoFs) = bcVector(edgeDoFs) + edgeVector;
    end
end

end