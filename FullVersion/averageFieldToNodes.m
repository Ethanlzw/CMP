function nodalField = averageFieldToNodes(model, results, fieldName)
% Average element field values to nodes

if ~isfield(results, fieldName)
    error('Field "%s" is not found in results.', fieldName);
end

elementFieldData = results.(fieldName);
nElements = size(elementFieldData,1);
nComponents = size(elementFieldData,2);

nNodes = size(model.nodesInfo,1);

nodalField = zeros(nNodes, nComponents);
nodeContributionCount = zeros(nNodes, 1);

elementTypeTable = model.elesType;
elementSetNames = model.elesSetName;
nElementSets = size(elementTypeTable,1);

iResult = 0;
for iSet = 1:nElementSets
    nElementNodes = elementTypeTable{iSet,4};
    elementSet = model.elesInfo.(elementSetNames{iSet});

    for iElement = 1:size(elementSet,1)
        iResult = iResult + 1;
        if iResult > nElements
            error('results.%s size is inconsistent with model elements.', fieldName);
        end

        elementNodeIDs = elementSet(iElement,2:1+nElementNodes);
        currentValue = elementFieldData(iResult,:);

        for a = 1:nElementNodes
            nodeID = elementNodeIDs(a);
            nodalField(nodeID,:) = nodalField(nodeID,:) + currentValue;
            nodeContributionCount(nodeID) = nodeContributionCount(nodeID) + 1;
        end
    end
end

if iResult ~= nElements
    error('results.%s size is inconsistent with model elements.', fieldName);
end

for iNode = 1:nNodes
    if nodeContributionCount(iNode) > 0
        nodalField(iNode,:) = nodalField(iNode,:) / nodeContributionCount(iNode);
    end
end

end