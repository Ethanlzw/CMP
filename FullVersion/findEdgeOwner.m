function edgeInfo = findEdgeOwner(model, nodeA, nodeB)
% Find one element edge containing boundary edge (nodeA,nodeB)

elementTypeTable = model.elesType;
elementSetNames = model.elesSetName;
nElementSets = size(elementTypeTable,1);

edgeInfo = struct();
edgeInfo.found = false;

for iSet = 1:nElementSets
    elementName = elementTypeTable{iSet,2};
    nElementNodes = elementTypeTable{iSet,4};
    elementSet = model.elesInfo.(elementSetNames{iSet});
    nElementsThisSet = size(elementSet,1);

    for iElement = 1:nElementsThisSet
        elementID = elementSet(iElement,1);
        elementNodeIDs = elementSet(iElement,2:1+nElementNodes);

        switch nElementNodes
            case 3
                edgeList = [elementNodeIDs([1 2]);
                            elementNodeIDs([2 3]);
                            elementNodeIDs([3 1])];
            case 4
                edgeList = [elementNodeIDs([1 2]);
                            elementNodeIDs([2 3]);
                            elementNodeIDs([3 4]);
                            elementNodeIDs([4 1])];
            otherwise
                error('Unsupported nEleNodes = %d', nElementNodes);
        end

        for k = 1:size(edgeList,1)
            edgeNode1 = edgeList(k,1);
            edgeNode2 = edgeList(k,2);

            if (edgeNode1 == nodeA && edgeNode2 == nodeB) || (edgeNode1 == nodeB && edgeNode2 == nodeA)
                edgeInfo.found = true;
                edgeInfo.setI = iSet;
                edgeInfo.eleID = elementID;
                edgeInfo.eleName = elementName;
                edgeInfo.edgeNodes = [edgeNode1 edgeNode2];
                edgeInfo.allNodeIDs = elementNodeIDs;
                return;
            end
        end
    end
end

end