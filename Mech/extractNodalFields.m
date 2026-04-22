function fields = extractNodalFields(model, u)
% extractNodalFields  Split global solution vector into per-field nodal arrays.
%
%   Output fields.(label) is a [nNodes × 1] vector for each DOF label.

dofInfo  = getDofInfo(model);
nNodes   = size(model.nodesInfo, 1);
nDof     = dofInfo.nDofPerNode;

% Reshape u into [nDof × nNodes] and extract each row
U = reshape(u, nDof, nNodes);   % row i = DOF i for all nodes

fields = struct();
for i = 1:nDof
    fields.(dofInfo.labels{i}) = U(i, :)';
end

end
