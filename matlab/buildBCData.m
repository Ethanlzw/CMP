function bcData = buildBCData(model)
% Convert raw BC definitions into unified data structure

dofInfo = getDofInfo(model);
nodesNum = size(model.nodesInfo, 1);
totalDofs = nodesNum * dofInfo.nDofPerNode;

bcData = struct();
bcData.totalDofs = totalDofs;
bcData.essential = struct();
bcData.natural   = struct();

%% Essential BC
uPrescribed = nan(totalDofs, 1);

if isfield(model, 'bcDef') && isfield(model.bcDef, 'essential')
    ess = model.bcDef.essential;

    if ~isempty(ess) && isfield(ess, 'nodes')
        for k = 1:numel(ess.nodes)
            nodeId = ess.nodes(k);
            fieldName = ess.fields{k};
            value = ess.values(k);

            localId = dofInfo.map.(fieldName);
            globalId = (nodeId - 1) * dofInfo.nDofPerNode + localId;
            uPrescribed(globalId) = value;
        end
    end
end

bcData.essential.uPrescribed = uPrescribed;
bcData.essential.knownDofs   = find(~isnan(uPrescribed));
bcData.essential.unknownDofs = find(isnan(uPrescribed));

%% Natural BC
if isfield(model, 'bcDef') && isfield(model.bcDef, 'natural')
    nat = model.bcDef.natural;
    bcData.natural = nat;
else
    bcData.natural.edges   = [];
    bcData.natural.type    = {};
    bcData.natural.values  = [];
    bcData.natural.setName = {};
end

end