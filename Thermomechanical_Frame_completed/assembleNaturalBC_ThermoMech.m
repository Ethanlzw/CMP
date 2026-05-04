function [bcMatrix, bcVector] = assembleNaturalBC_ThermoMech(model)
% Assemble thermomechanical natural BC contributions from model.natBC
%
% model.natBC rows:
%   { [node1,node2,...], [tx, ty, M, S] }

nNodeDoF = model.elesType{1,3};
if nNodeDoF ~= 3
    error('assembleNaturalBC_ThermoMech is only for models with 3 DoFs per node.');
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

    if numel(bcVals) ~= 4
        error('Thermomechanical natural BC values must be [tx, ty, M, S].');
    end

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

        [edgeMatrix, edgeVector] = Edge2_ThermoMech(nodeCoords, bcVals, thickness);

        %% Complete the following lines
        edgeDoFs = [(nodeA-1)*nNodeDoF+1;
                    (nodeA-1)*nNodeDoF+2;
                    (nodeA-1)*nNodeDoF+3;
                    (nodeB-1)*nNodeDoF+1;
                    (nodeB-1)*nNodeDoF+2;
                    (nodeB-1)*nNodeDoF+3];

        bcMatrix(edgeDoFs, edgeDoFs) = bcMatrix(edgeDoFs, edgeDoFs) + edgeMatrix;
        bcVector(edgeDoFs) = bcVector(edgeDoFs) + edgeVector;
    end
end

end