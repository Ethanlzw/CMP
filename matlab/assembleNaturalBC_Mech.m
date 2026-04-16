function bcVector = assembleNaturalBC_Mech(model)
% Assemble mechanical natural BC contributions

dofInfo = getDofInfo(model);
nodesNum = size(model.nodesInfo, 1);
bcVector = zeros(nodesNum * dofInfo.nDofPerNode, 1);

if ~isfield(model, 'bcData') || ~isfield(model.bcData, 'natural')
    return;
end

nat = model.bcData.natural;

if isempty(nat) || isempty(nat.edges)
    return;
end

for i = 1:size(nat.edges, 1)
    edgeNodes = nat.edges(i, :);
    edgeNodesInfo = model.nodesInfo(edgeNodes, :);

    bc = struct();
    bc.type = nat.type{i};
    bc.value = nat.values(i, :).';

    state = struct();

    [~, edgeVector] = Edge2_MechTraction(edgeNodesInfo, bc, model.problem, model.loads, state);
    edgeDof = getElementDof(edgeNodes, dofInfo);

    bcVector(edgeDof) = bcVector(edgeDof) + edgeVector;
end

end