function [u, reaction] = solveSystem(K, bcData, F)
% solveSystem  Apply essential BCs via partitioning and solve K*u = F.
%
%   Partitions the system into known (k) and unknown (u) DOFs:
%       K_uu * u_u = F_u - K_uk * u_k
%
%   Outputs:
%     u        – full displacement vector
%     reaction – nodal reaction vector  R = K*u - F

knownDofs   = bcData.essential.knownDofs;
unknownDofs = bcData.essential.unknownDofs;
uPrescribed = bcData.essential.uPrescribed;

u = zeros(size(F));
u(knownDofs) = uPrescribed(knownDofs);

% Solve for free DOFs
Kuu = K(unknownDofs, unknownDofs);
Kuk = K(unknownDofs, knownDofs);
u(unknownDofs) = Kuu \ (F(unknownDofs) - Kuk * u(knownDofs));

reaction = K * u - F;

end
