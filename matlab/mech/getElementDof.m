function eleDof = getElementDof(eleNodes, dofInfo)
% getElementDof  Return global DOF indices for one element (vectorised).
%
%   Layout: [node1_dof1, node1_dof2, ..., node2_dof1, ...]

nDof   = dofInfo.nDofPerNode;
base   = (eleNodes(:) - 1) * nDof;            % [nEleNodes × 1]
offsets = 1:nDof;                              % [1 × nDof]
eleDof = (base + offsets)';                    % [nDof × nEleNodes], then linearise
eleDof = eleDof(:)';                           % [1 × nEleNodes*nDof]

end
