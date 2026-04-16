function fields = extractNodalFields(model, u)
% Extract nodal field arrays from global solution vector

dofInfo = getDofInfo(model);
nodesNum = size(model.nodesInfo, 1);

fields = struct();

for i = 1:numel(dofInfo.labels)
    label = dofInfo.labels{i};
    fields.(label) = zeros(nodesNum, 1);
end

for nodeId = 1:nodesNum
    base = (nodeId - 1) * dofInfo.nDofPerNode;
    for i = 1:numel(dofInfo.labels)
        label = dofInfo.labels{i};
        fields.(label)(nodeId) = u(base + i);
    end
end

end