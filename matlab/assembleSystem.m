function [globalMatrix, globalVector] = assembleSystem(model)
% Assemble global matrix and global vector for mixed element sets

nodesInfo   = model.nodesInfo;
elesType    = model.elesType;
elesSetName = model.elesSetName;

dofInfo = getDofInfo(model);

nodesNum   = size(nodesInfo, 1);
elesSetNum = size(elesType, 1);

nNodeDoF = elesType{1,3};
for setI = 2:elesSetNum
    if elesType{setI,3} ~= nNodeDoF
        error('All element sets in one model must have the same node DoF count.');
    end
end

globalMatrix = zeros(nodesNum * nNodeDoF, nodesNum * nNodeDoF);
globalVector = zeros(nodesNum * nNodeDoF, 1);

for setI = 1:elesSetNum
    eleName = elesType{setI,2};
    eleData = getElementData(eleName, model);
    eleFunc = eleData.func;

    elesSet = model.elesInfo.(elesSetName{setI});
    if isempty(elesSet)
        continue;
    end

    nSetEles = size(elesSet, 1);

    for eleI = 1:nSetEles
        eleNodes = elesSet(eleI, :);
        eleNodesInfo = nodesInfo(eleNodes, :);
        eleDof = getElementDof(eleNodes, dofInfo);

        material = model.material;
        problem  = model.problem;
        loads    = model.loads;
        state    = struct();

        [eleMatrix, eleVector] = eleFunc(eleNodesInfo, material, problem, loads, state);

        globalMatrix(eleDof, eleDof) = globalMatrix(eleDof, eleDof) + eleMatrix;
        globalVector(eleDof) = globalVector(eleDof) + eleVector;
    end
end

end