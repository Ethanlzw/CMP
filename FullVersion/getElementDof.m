function elementDofs = getElementDof(elementNodeIDs, dofInfo)
% Return global DoF indices of one element

nNodeDofs = dofInfo.nDoFPerNode;
nElementNodes = numel(elementNodeIDs);

elementDofs = zeros(nElementNodes*nNodeDofs,1);

for a = 1:nElementNodes
    for iDof = 1:nNodeDofs
        elementDofs((a-1)*nNodeDofs + iDof) = (elementNodeIDs(a)-1)*nNodeDofs + iDof;
    end
end

end