function [reaction, solution] = solveSystem(globalMatrix, globalVector, bcData)
% Solve linear FE system with prescribed DoFs

nTotalDofs = length(globalVector);
allDofs = 1:nTotalDofs;

prescribedDofs = bcData.prescribedDoFs;
prescribedValues = bcData.prescribedValues;
activeDofs = setdiff(allDofs, prescribedDofs);

solution = zeros(nTotalDofs,1);
solution(prescribedDofs) = prescribedValues(prescribedDofs);

solution(activeDofs) = ...
    globalMatrix(activeDofs,activeDofs) \ ...
   (globalVector(activeDofs) - globalMatrix(activeDofs,prescribedDofs) * solution(prescribedDofs));

reaction = globalMatrix * solution - globalVector;

end