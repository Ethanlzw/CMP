function [K, F] = assembleSystem(model)
% assembleSystem  Assemble global stiffness matrix K and internal load vector F.
%
%   Iterates over all element sets, calls element routines, and accumulates
%   contributions into the global system via direct stiffness assembly.

nodesInfo = model.nodesInfo;
elesSets  = model.elesSets;       % {setField, elementName; ...}
dofInfo   = getDofInfo(model);

nNodes  = size(nodesInfo, 1);
nTotalDof = nNodes * dofInfo.nDofPerNode;

K = zeros(nTotalDof, nTotalDof);
F = zeros(nTotalDof, 1);

for s = 1:size(elesSets, 1)
    setField = elesSets{s, 1};
    eleName  = elesSets{s, 2};

    elesSet = model.elesInfo.(setField);
    if isempty(elesSet), continue; end

    eleData = getElementData(eleName, model);
    eleFunc = eleData.func;

    for e = 1:size(elesSet, 1)
        eleNodes     = elesSet(e, :);
        eleNodesInfo = nodesInfo(eleNodes, :);
        eleDof       = getElementDof(eleNodes, dofInfo);

        [Ke, Fe] = eleFunc(eleNodesInfo, model.material, model.problem, model.loads, struct());

        K(eleDof, eleDof) = K(eleDof, eleDof) + Ke;
        F(eleDof)         = F(eleDof)         + Fe;
    end
end

end
