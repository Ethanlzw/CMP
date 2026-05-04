function results = postprocessModel(model, solution)
% Postprocess model with mixed element sets
%
% Output:
%   results(e).id
%   results(e).name
%   results(e).set
%   results(e).nodeIDs
%   results(e).value

nodesInfo = model.nodesInfo;
elesType = model.elesType;
elesSetName = model.elesSetName;

elesSetNum = size(elesType,1);
nNodeDoF = elesType{1,3};

if isfield(model,'problem')
    problem = model.problem;
else
    problem = struct();
end

results = struct([]);
resCount = 0;

for setI = 1:elesSetNum
    eleName = elesType{setI,2};
    nEleNodes = elesType{setI,4};

    elesSet = model.elesInfo.(elesSetName{setI});
    nSetEles = size(elesSet,1);

    postFunc = str2func(['post_' eleName]);

    for e = 1:nSetEles
        resCount = resCount + 1;

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

        eleDoFs = getElementDoFs(eleNodeIDs, nNodeDoF);
        eleSol = solution(eleDoFs);

        results(resCount).id = eleID;
        results(resCount).name = eleName;
        results(resCount).set = setI;
        results(resCount).nodeIDs = eleNodeIDs;
        results(resCount).value = postFunc(eleNodesInfo, material, problem, eleSol);
    end
end

end