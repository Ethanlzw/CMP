function edgeInfo = findEdgeOwner(model, nodeA, nodeB)
% Find one element edge containing boundary edge (nodeA,nodeB)

elesType = model.elesType;
elesSetName = model.elesSetName;
elesSetNum = size(elesType,1);

edgeInfo = struct();
edgeInfo.found = false;

for setI = 1:elesSetNum
    eleName = elesType{setI,2};
    nEleNodes = elesType{setI,4};
    elesSet = model.elesInfo.(elesSetName{setI});
    nSetEles = size(elesSet,1);

    for e = 1:nSetEles
        eleID = elesSet(e,1);
        eleNodeIDs = elesSet(e,2:1+nEleNodes);

        % build local edges
        switch nEleNodes
            case 3
                edgeList = [eleNodeIDs([1 2]);
                            eleNodeIDs([2 3]);
                            eleNodeIDs([3 1])];
            case 4
                edgeList = [eleNodeIDs([1 2]);
                            eleNodeIDs([2 3]);
                            eleNodeIDs([3 4]);
                            eleNodeIDs([4 1])];
            otherwise
                error('Unsupported nEleNodes = %d', nEleNodes);
        end

        for k = 1:size(edgeList,1)
            e1 = edgeList(k,1);
            e2 = edgeList(k,2);

            if (e1 == nodeA && e2 == nodeB) || (e1 == nodeB && e2 == nodeA)
                edgeInfo.found = true;
                edgeInfo.setI = setI;
                edgeInfo.eleID = eleID;
                edgeInfo.eleName = eleName;
                edgeInfo.edgeNodes = [e1 e2];
                edgeInfo.allNodeIDs = eleNodeIDs;
                return;
            end
        end
    end
end

end