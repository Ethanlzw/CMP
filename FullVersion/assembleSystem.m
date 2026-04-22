function [globalMatrix, globalVector] = assembleSystem(model)
% Assemble domain contributions from all element sets

nodesInfo = model.nodesInfo;
elementTypeTable = model.elesType;
elementSetNames = model.elesSetName;

dofInfo = getDofInfo(model);

nNodes = size(nodesInfo,1);
nElementSets = size(elementTypeTable,1);
nNodeDofs = dofInfo.nDoFPerNode;

for iSet = 1:nElementSets
    if elementTypeTable{iSet,3} ~= nNodeDofs
        error('All element sets in one model must have the same node DoF count.');
    end
end

globalMatrix = zeros(nNodes*nNodeDofs, nNodes*nNodeDofs);
globalVector = zeros(nNodes*nNodeDofs, 1);

for iSet = 1:nElementSets
    elementName = elementTypeTable{iSet,2};
    elementInfo = getElementData(elementName, model);
    nElementNodes = elementInfo.nEleNodes;

    elementSet = model.elesInfo.(elementSetNames{iSet});
    if isempty(elementSet)
        continue;
    end

    nElementsThisSet = size(elementSet,1);

    for iElement = 1:nElementsThisSet
        elementID = elementSet(iElement,1);
        elementNodeIDs = elementSet(iElement,2:1+nElementNodes);

        elementNodesInfo = zeros(nElementNodes, size(nodesInfo,2));
        for a = 1:nElementNodes
            nodeMask = nodesInfo(:,1) == elementNodeIDs(a);

            if nnz(nodeMask) == 0
                error('Node ID %d of element %d is not found in model.nodesInfo.', elementNodeIDs(a), elementID);
            elseif nnz(nodeMask) > 1
                error('Node ID %d of element %d appears more than once in model.nodesInfo.', elementNodeIDs(a), elementID);
            end

            elementNodesInfo(a,:) = nodesInfo(nodeMask,:);
        end

        elementDofs = getElementDof(elementNodeIDs, dofInfo);

        sectionMask = model.section(:,1) == elementID;
        if nnz(sectionMask) == 0
            error('Element %d is not found in model.section.', elementID);
        elseif nnz(sectionMask) > 1
            error('Element %d appears more than once in model.section.', elementID);
        end
        materialData = model.section(sectionMask,2:end);

        if isfield(model,'problem')
            problemData = model.problem;
        else
            problemData = struct();
        end

        if isfield(model,'state')
            stateData = model.state;
        else
            stateData = struct();
        end

        loadData = struct();
        if isfield(model,'eleLoadData')
            if numel(model.eleLoadData) >= elementID && ~isempty(model.eleLoadData{elementID})
                loadData = model.eleLoadData{elementID};
            end
        end

        [elementMatrix, elementVector] = elementInfo.func(elementNodesInfo, materialData, problemData, loadData, stateData);

        globalMatrix(elementDofs, elementDofs) = globalMatrix(elementDofs, elementDofs) + elementMatrix;
        globalVector(elementDofs) = globalVector(elementDofs) + elementVector;
    end
end

if isfield(model,'nodalVector') && ~isempty(model.nodalVector)
    globalVector = globalVector + model.nodalVector;
end

end