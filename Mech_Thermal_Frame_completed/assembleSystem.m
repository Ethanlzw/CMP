function [globalMatrix, globalVector] = assembleSystem(model)
% Assemble global matrix and global vector for mixed element sets
% 主逻辑:
%   1) 遍历每个单元集
%   2) 读取单元节点 / 材料 / 载荷
%   3) 调用对应单元例程得到局部矩阵和向量
%   4) 通过自由度映射装配到整体系统
%
% 这里仅处理域内单元项.
% 本质边界条件不在这里施加.
% 自然边界条件由 buildSystem 中的 assembleNaturalBC_xxx 处理.
nodesInfo = model.nodesInfo;
elesType = model.elesType;
elesSetName = model.elesSetName;

nodesNum = size(nodesInfo,1);
elesSetNum = size(elesType,1);

%% 所有单元集必须共享同样的每节点自由度数
nNodeDoF = elesType{1,3};
for setI = 2:elesSetNum
    if elesType{setI,3} ~= nNodeDoF
        error('All element sets in one model must have the same node DoF count.');
    end
end

globalMatrix = zeros(nodesNum*nNodeDoF, nodesNum*nNodeDoF);
globalVector = zeros(nodesNum*nNodeDoF, 1);

for setI = 1:elesSetNum
    eleName = elesType{setI,2};
    eleData = getElementData(eleName, model);
    nEleNodes = eleData.nEleNodes;

    elesSet = model.elesInfo.(elesSetName{setI});
    if isempty(elesSet)
        continue;
    end

    nSetEles = size(elesSet,1);

    for e = 1:nSetEles
        eleID = elesSet(e,1);
        eleNodeIDs = elesSet(e,2:1+nEleNodes);
        
        % 取出当前单元的节点坐标信息
        eleNodesInfo = zeros(nEleNodes, size(nodesInfo,2));
        for a = 1:nEleNodes
            nodeMask = nodesInfo(:,1) == eleNodeIDs(a);
            eleNodesInfo(a,:) = nodesInfo(nodeMask,:);
        end
        
        % 按 element ID 读取材料 / 截面数据
        matMask = model.section(:,1) == eleID;
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
        
        % 默认单元载荷为空结构体
        loads = struct();
        if isfield(model,'eleLoadData')
            if numel(model.eleLoadData) >= eleID && ~isempty(model.eleLoadData{eleID})
                loads = model.eleLoadData{eleID};
            end
        end
        
        %% CORE LOGIC: 单元函数调用, DoF 映射, 整体矩阵与整体向量装配
        % 调用单元
        [eleMatrix, eleVector] = eleData.func(eleNodesInfo, material, problem, loads, state);

        %% ----- Add three lines here to assemble the global matrix and vector
        %(without boundary conditions)

        % 先得到当前单元的全局自由度列表
        % 当前单元局部自由度到全局自由度的映射
        eleDoFs = getElementDoFs(eleNodeIDs, nNodeDoF);

        % 装配局部刚度矩阵与局部载荷向量
        globalMatrix(eleDoFs, eleDoFs) = globalMatrix(eleDoFs, eleDoFs) + eleMatrix;
        globalVector(eleDoFs) = globalVector(eleDoFs) + eleVector;

    end
end

% 若模型中有直接施加的全局节点力, 在最后统一叠加
if isfield(model,'nodalVector') && ~isempty(model.nodalVector)
    globalVector = globalVector + model.nodalVector;
end

end