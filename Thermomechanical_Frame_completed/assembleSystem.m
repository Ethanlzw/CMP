function [globalMatrix, globalVector] = assembleSystem(model)
% Assemble global matrix and global vector for mixed element sets

nodesInfo = model.nodesInfo;
elesType = model.elesType;
elesSetName = model.elesSetName;

nodesNum = size(nodesInfo,1);
elesSetNum = size(elesType,1);

% All element sets in one model must have the same number of DoFs per node
nNodeDoF = elesType{1,3};
for setI = 2:elesSetNum
    if elesType{setI,3} ~= nNodeDoF
        error('All element sets in one model must have the same node DoF count.');
    end
end

% Initialize global system
globalMatrix = zeros(nodesNum*nNodeDoF, nodesNum*nNodeDoF);
globalVector = zeros(nodesNum*nNodeDoF, 1);

% Loop over element sets
for setI = 1:elesSetNum
    % Read this element set first
    elesSet = model.elesInfo.(elesSetName{setI});

    % Empty element set contributes nothing to the global system
    if isempty(elesSet)
        continue;
    end

    % Only read element data when this set really contains elements
    eleName = elesType{setI,2};
    eleData = getElementData(eleName, model);
    nEleNodes = eleData.nEleNodes;

    nSetEles = size(elesSet,1);

    % Loop over elements in this set
    for e = 1:nSetEles
        eleID = elesSet(e,1);
        eleNodeIDs = elesSet(e,2:1+nEleNodes);

        % Collect nodal information of the current element
        eleNodesInfo = zeros(nEleNodes, size(nodesInfo,2));
        for a = 1:nEleNodes
            nodeMask = nodesInfo(:,1) == eleNodeIDs(a);
            eleNodesInfo(a,:) = nodesInfo(nodeMask,:);
        end

        % Read material / section data of the current element
        matMask = model.section(:,1) == eleID;
        material = model.section(matMask,2:end);
        
        % Read optional problem data
        if isfield(model,'problem')
            problem = model.problem;
        else
            problem = struct();
        end

        % Read optional state data
        if isfield(model,'state')
            state = model.state;
        else
            state = struct();
        end

        % Read optional element load data
        loads = struct();
        if isfield(model,'eleLoadData')
            if numel(model.eleLoadData) >= eleID && ~isempty(model.eleLoadData{eleID})
                loads = model.eleLoadData{eleID};
            end
        end

        % Compute element matrix and element vector
        [eleMatrix, eleVector] = eleData.func(eleNodesInfo, material, problem, loads, state);

        %----- Add three lines here to assemble the global matrix and vector
        %(without boundary conditions)
        % Assemble into global matrix and global vector
        eleDoFs = getElementDoFs(eleNodeIDs, nNodeDoF);
        globalMatrix(eleDoFs, eleDoFs) = globalMatrix(eleDoFs, eleDoFs) + eleMatrix;
        globalVector(eleDoFs) = globalVector(eleDoFs) + eleVector;

    end
end

% Add nodal force / source vector if provided
if isfield(model,'nodalVector') && ~isempty(model.nodalVector)
    globalVector = globalVector + model.nodalVector;
end

end