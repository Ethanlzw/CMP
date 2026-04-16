function result = postprocessModel(model, u)
% Unified post-processing entry

dofInfo = getDofInfo(model);
result = struct();

result.u = u;
result.fields = extractNodalFields(model, u);
result.element = struct();

for setI = 1:numel(model.elesSetName)
    setName = model.elesSetName{setI};
    eleName = model.elesType{setI,2};
    elesSet = model.elesInfo.(setName);

    if isempty(elesSet)
        continue;
    end

    postFuncName = ['post_' eleName];
    postFunc = str2func(postFuncName);

    outCell = cell(size(elesSet,1), 1);

    for eleI = 1:size(elesSet,1)
        eleNodes = elesSet(eleI,:);
        eleNodesInfo = model.nodesInfo(eleNodes,:);
        eleDof = getElementDof(eleNodes, dofInfo);
        eleU = u(eleDof);

        state = struct();
        outCell{eleI} = postFunc(eleNodesInfo, eleU, model.material, model.problem, state);
    end

    result.element.(setName) = outCell;
end

end