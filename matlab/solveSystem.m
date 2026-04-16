function [u, reaction] = solveSystem(K, bcData, F)
% Solve linear system with essential BC

uPrescribed = bcData.essential.uPrescribed;
knownDofs   = bcData.essential.knownDofs;
unknownDofs = bcData.essential.unknownDofs;

u = nan(size(F));
u(knownDofs) = uPrescribed(knownDofs);

Kuu = K(unknownDofs, unknownDofs);
Kuk = K(unknownDofs, knownDofs);

Fu = F(unknownDofs);
uk = u(knownDofs);

u(unknownDofs) = Kuu \ (Fu - Kuk * uk);

reaction = K * u - F;

end