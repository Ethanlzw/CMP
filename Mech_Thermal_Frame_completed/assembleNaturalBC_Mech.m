function bcVector = assembleNaturalBC_Mech(model)
% Assemble mechanical natural BC contributions from model.natBC
%
% model.natBC rows:
%   { [node1,node2,...], [tx,ty] }

nNodeDoF = model.elesType{1,3};
if nNodeDoF ~= 2
    error('assembleNaturalBC_Mech is only for mechanical models with 2 DoFs per node.');
end

numNodes = size(model.nodesInfo,1);
bcVector = zeros(numNodes*nNodeDoF,1);

if ~isfield(model,'natBC') || isempty(model.natBC)
    return;
end

if isfield(model,'problem') && isfield(model.problem,'thickness')
    thickness = model.problem.thickness;
else
    thickness = 1.0;
end

%% CORE LOGIC: 边界链拆边、节点坐标提取、Edge2 装配到全局向量
for ibc = 1:size(model.natBC,1)
    bcNodes = model.natBC{ibc,1};
    bcVals  = model.natBC{ibc,2};

    if numel(bcVals) ~= 2
        error('Mechanical natural BC values must be [tx, ty].');
    end
    
     % 一条边界链被拆成若干线单元, 每个线单元由相邻两节点定义
    for i = 1:(numel(bcNodes)-1)
        nodeA = bcNodes(i);
        nodeB = bcNodes(i+1);

        % 先确认这条边确实属于网格中的某个单元边.
        edgeInfo = findEdgeOwner(model, nodeA, nodeB);
        if ~edgeInfo.found
            error('Boundary edge (%d,%d) was not found in any element.', nodeA, nodeB);
        end

%% need write all these (find out everything)

        idxA = find(model.nodesInfo(:,1) == nodeA, 1);
        idxB = find(model.nodesInfo(:,1) == nodeB, 1);

        nodeCoords = [model.nodesInfo(idxA, 2:3);
                      model.nodesInfo(idxB, 2:3)];

        edgeVector = Edge2_MechTraction(nodeCoords, bcVals, thickness);

        % edgeDoFs = [; ; ; ];
     
        edgeDoFs = [(nodeA-1)*nNodeDoF+1;
                    (nodeA-1)*nNodeDoF+2;
                    (nodeB-1)*nNodeDoF+1;
                    (nodeB-1)*nNodeDoF+2];

        % 装配到全局载荷向量
        bcVector(edgeDoFs) = bcVector(edgeDoFs) + edgeVector;
    end
end

end