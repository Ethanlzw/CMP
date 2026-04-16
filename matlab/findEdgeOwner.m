function owner = findEdgeOwner(model, edgeNodes)
% Find which element set / element owns a given edge

owner = struct();
owner.found = false;
owner.setName = '';
owner.eleId = [];
owner.localFaceId = [];

edgeNodes = sort(edgeNodes(:)).';

for setI = 1:numel(model.elesSetName)
    setName = model.elesSetName{setI};
    eleName = model.elesType{setI,2};
    elesSet = model.elesInfo.(setName);

    if isempty(elesSet)
        continue;
    end

    faces = getElementFaces2D(eleName);

    for eleI = 1:size(elesSet,1)
        conn = elesSet(eleI,:);
        for fI = 1:size(faces,1)
            testEdge = sort(conn(faces(fI,:)));
            if isequal(testEdge, edgeNodes)
                owner.found = true;
                owner.setName = setName;
                owner.eleId = eleI;
                owner.localFaceId = fI;
                return;
            end
        end
    end
end

end