function bcData = buildBCData(model)
% Build prescribed DoF data from model.essBC
%
% model.essBC rows:
%   { [node1,node2,...], [v1,v2,...] }
%
% Any NaN entry means that DoF is not prescribed.

if ~isfield(model,'fieldNames')
    error('model.fieldNames must be defined.');
end

dofInfo = getDofInfo(model);

nNodes = size(model.nodesInfo,1);
nNodeDofs = dofInfo.nDoFPerNode;
nTotalDofs = nNodes * nNodeDofs;

prescribedMask = false(nTotalDofs,1);
prescribedValues = zeros(nTotalDofs,1);

if ~isfield(model,'essBC') || isempty(model.essBC)
    bcData = struct();
    bcData.prescribedDoFs = [];
    bcData.prescribedValues = prescribedValues;
    bcData.activeDoFs = (1:nTotalDofs).';
    return;
end

valueTol = 1e-12;

for iBC = 1:size(model.essBC,1)
    bcNodeIDs = model.essBC{iBC,1};
    bcValuesPerNode = model.essBC{iBC,2};

    if numel(bcValuesPerNode) ~= nNodeDofs
        error('Essential BC values must have length %d.', nNodeDofs);
    end

    for iNode = 1:numel(bcNodeIDs)
        nodeID = bcNodeIDs(iNode);

        if nodeID < 1 || nodeID > nNodes
            error('Node ID %d in model.essBC is out of range.', nodeID);
        end

        for iDof = 1:nNodeDofs
            if ~isnan(bcValuesPerNode(iDof))
                globalDof = (nodeID - 1) * nNodeDofs + iDof;

                if prescribedMask(globalDof)
                    if abs(prescribedValues(globalDof) - bcValuesPerNode(iDof)) > valueTol
                        error('Conflicting essential BC at global DoF %d.', globalDof);
                    end
                else
                    prescribedMask(globalDof) = true;
                    prescribedValues(globalDof) = bcValuesPerNode(iDof);
                end
            end
        end
    end
end

bcData = struct();
bcData.prescribedDoFs = find(prescribedMask);
bcData.prescribedValues = prescribedValues;
bcData.activeDoFs = find(~prescribedMask);

end