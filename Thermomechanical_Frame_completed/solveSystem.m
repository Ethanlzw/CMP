function [reaction, solution] = solveSystem(globalMatrix, globalVector, bcData)
% Solve linear FE system with prescribed DoFs

totalDoFs = length(globalVector);
allDoFs = 1:totalDoFs;

prescribedDoFs = bcData.prescribedDoFs;
prescribedValues = bcData.prescribedValues;

activeDoFs = setdiff(allDoFs, prescribedDoFs);

solution = zeros(totalDoFs,1);
solution(prescribedDoFs) = prescribedValues(prescribedDoFs);

%% Complete the following session

Kaa = globalMatrix(activeDoFs, activeDoFs);
Kap = globalMatrix(activeDoFs, prescribedDoFs);
Fa = globalVector(activeDoFs);
up = solution(prescribedDoFs);

solution(activeDoFs) = Kaa \ (Fa - Kap * up);

reaction = globalMatrix * solution - globalVector;

end