function bcData = buildBCData(model)
% buildBCData  Convert model.bcDef into solver-ready BC data.
%
%   bcData.essential.uPrescribed  – full-length vector, NaN for free DOFs
%   bcData.essential.knownDofs    – indices of prescribed DOFs
%   bcData.essential.unknownDofs  – indices of free DOFs
%   bcData.forceVector            – external nodal force vector

dofInfo   = getDofInfo(model);
nNodes    = size(model.nodesInfo, 1);
totalDofs = nNodes * dofInfo.nDofPerNode;

% ── Essential (Dirichlet) BCs ─────────────────────────────────────────────
uPrescribed = nan(totalDofs, 1);

if isfield(model, 'bcDef') && isfield(model.bcDef, 'essential')
    ess     = model.bcDef.essential;
    nodeIds = ess.nodes(:);
    values  = ess.values(:);

    % Vectorised global DOF computation
    localDofs = cellfun(@(f) dofInfo.map.(f), ess.fields(:));
    globalDofs = (nodeIds - 1) * dofInfo.nDofPerNode + localDofs;

    uPrescribed(globalDofs) = values;
end

bcData.essential.uPrescribed = uPrescribed;
bcData.essential.knownDofs   = find(~isnan(uPrescribed));
bcData.essential.unknownDofs = find( isnan(uPrescribed));

% ── Natural (Neumann) BCs ─────────────────────────────────────────────────
if isfield(model, 'bcDef') && isfield(model.bcDef, 'forceVector')
    bcData.forceVector = model.bcDef.forceVector(:);
else
    bcData.forceVector = zeros(totalDofs, 1);
end

end
