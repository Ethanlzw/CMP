function nodal = extractNodalFields(model, solution)
% Extract nodal fields from global solution vector

dofInfo = getDofInfo(model);

nNodes = size(model.nodesInfo,1);
nNodeDofs = dofInfo.nDoFPerNode;

if numel(solution) ~= nNodes * nNodeDofs
    error('Solution length is inconsistent with model size.');
end

nodal = struct();

for iDof = 1:nNodeDofs
    fieldLabel = dofInfo.labels{iDof};
    nodal.(fieldLabel) = solution(iDof:nNodeDofs:end);
end

end