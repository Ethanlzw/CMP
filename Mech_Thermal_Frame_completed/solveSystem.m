function [reaction, solution] = solveSystem(globalMatrix, globalVector, bcData)
% Solve linear FE system with prescribed DoFs
% 方程形式:
%   K * u = F
%
% 处理方法:
%   将自由度分成
%     - prescribedDoFs : 已知位移 / 温度等本质边界
%     - activeDoFs     : 需要求解的自由度
%
% 分块后有:
%   Kaa * ua + Kap * up = Fa
%   -> ua = Kaa \ (Fa - Kap * up)
%
% 输出:
%   reaction : 全部自由度上的反力 / 平衡残量
%   solution : 全局解向量

totalDoFs = length(globalVector);
allDoFs = 1:totalDoFs;

prescribedDoFs = bcData.prescribedDoFs;
prescribedValues = bcData.prescribedValues;

%% CORE LOGIC: 已知 / 未知自由度划分, 分块矩阵求解
% need understanding 
activeDoFs = setdiff(allDoFs, prescribedDoFs);

solution = zeros(totalDoFs,1);
% 写入已知自由度
solution(prescribedDoFs) = prescribedValues(prescribedDoFs);

Kaa = globalMatrix(activeDoFs, activeDoFs);
Kap = globalMatrix(activeDoFs, prescribedDoFs);
Fa = globalVector(activeDoFs);
up = solution(prescribedDoFs);

% 只在未知自由度子空间内求解
solution(activeDoFs) = Kaa \ (Fa - Kap * up);

% 反力统一由完整方程回代得到
reaction = globalMatrix * solution - globalVector;

end