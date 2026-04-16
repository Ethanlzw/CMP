function eleDof = getElementDof(eleNodes, dofInfo)
% Return global DoF indices for one element

nEleNodes = numel(eleNodes);
nDofPerNode = dofInfo.nDofPerNode;

eleDof = zeros(1, nEleNodes * nDofPerNode);

for i = 1:nEleNodes
    nodeId = eleNodes(i);
    base = (nodeId - 1) * nDofPerNode;
    idx = (i - 1) * nDofPerNode + (1:nDofPerNode);
    eleDof(idx) = base + (1:nDofPerNode);
end

end