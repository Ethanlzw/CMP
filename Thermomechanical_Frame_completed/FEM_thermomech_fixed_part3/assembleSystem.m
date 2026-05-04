function [globalMatrix, globalVector] = assembleSystem(model)
% Assemble global matrix and global vector for mixed element sets
% This version skips empty element sets before resolving element metadata.

nodesInfo = model.nodesInfo;
elesType = model.elesType;
elesSetName = model.elesSetName;

nodesNum = size(nodesInfo,1);
elesSetNum = size(elesType,1);

nNodeDoF = elesType{1,3};
for setI = 2:elesSetNum
    if elesType{setI,3} ~= nNodeDoF
        error('All element sets in one model must have the same node DoF count.');
    end
end

globalMatrix = zeros(nodesNum*nNodeDoF, nodesNum*nNodeDoF);
globalVector = zeros(nodesNum*nNodeDoF, 1);

for setI = 1:elesSetNum
    if ~isfield(model.elesInfo, elesSetName{setI})
        error('Element set "%s" is missing in model.elesInfo.', elesSetName{setI});
    end

    elesSet = model.elesInfo.(elesSetName{setI});
    if isempty(elesSet)
        continue;
    end

    eleName = elesType{setI,2};
    eleData = getElementData(eleName, model);
    nEleNodes = eleData.nEleNodes;
    nSetEles = size(elesSet,1);

    for e = 1:nSetEles
        eleID = elesSet(e,1);
        eleNodeIDs = elesSet(e,2:1+nEleNodes);

        eleNodesInfo = zeros(nEleNodes, size(nodesInfo,2));
        for a = 1:nEleNodes
            nodeMask = nodesInfo(:,1) == eleNodeIDs(a);
            if ~any(nodeMask)
                error('Node ID %d of element %d was not found.', eleNodeIDs(a), eleID);
            end
            eleNodesInfo(a,:) = nodesInfo(nodeMask,:);
        end

        matMask = model.section(:,1) == eleID;
        if ~any(matMask)
            error('Section row for element %d was not found.', eleID);
        end
        material = model.section(matMask,2:end);

        if isfield(model,'problem')
            problem = model.problem;
        else
            problem = struct();
        end

        if isfield(model,'state')
            state = model.state;
        else
            state = struct();
        end

        loads = struct();
        if isfield(model,'eleLoadData')
            if numel(model.eleLoadData) >= eleID && ~isempty(model.eleLoadData{eleID})
                loads = model.eleLoadData{eleID};
            end
        end

        [eleMatrix, eleVector] = eleData.func(eleNodesInfo, material, problem, loads, state);
        eleDoFs = getElementDoFs(eleNodeIDs, nNodeDoF);

        globalMatrix(eleDoFs, eleDoFs) = globalMatrix(eleDoFs, eleDoFs) + eleMatrix;
        globalVector(eleDoFs) = globalVector(eleDoFs) + eleVector;
    end
end

if isfield(model,'nodalVector') && ~isempty(model.nodalVector)
    globalVector = globalVector + model.nodalVector;
end

end
