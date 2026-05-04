function bcData = buildBCData(model)
% Build prescribed DoF data from model.essBC
%
% Supported format:
%   model.essBC = {
%       [node list], [value vector];
%       [node list], [value vector];
%       ...
%   }
%
% Examples:
%   Thermal (1 DoF/node):
%       model.essBC = {
%           [4,8,12], [-20];
%           [3,7,11], [0]
%       };
%
%   Mechanical (2 DoFs/node):
%       model.essBC = {
%           [1,4], [0,0];
%           [2],   [1e-3,0]
%       };
%
%   Thermomechanical (3 DoFs/node):
%       model.essBC = {
%           [1,4], [0,0,20];
%           [2],   [0,0,100]
%       };

nNodeDoF = model.elesType{1,3};
numNodes = size(model.nodesInfo,1);
totalDoFs = numNodes * nNodeDoF;

bcData = struct();
bcData.prescribedDoFs = [];
bcData.prescribedValues = zeros(totalDoFs,1);

isPrescribed = false(totalDoFs,1);

if ~isfield(model,'essBC') || isempty(model.essBC)
    return;
end

nBCGroups = size(model.essBC,1);

for ibc = 1:nBCGroups
    bcNodes = model.essBC{ibc,1};
    bcVals  = model.essBC{ibc,2};

    if numel(bcVals) ~= nNodeDoF
        error('Essential BC group %d has %d values, but nNodeDoF = %d.', ...
              ibc, numel(bcVals), nNodeDoF);
    end

    for i = 1:numel(bcNodes)
        nodeID = bcNodes(i);

        nodeMask = (model.nodesInfo(:,1) == nodeID);
        if ~any(nodeMask)
            error('Essential BC node %d was not found in nodesInfo.', nodeID);
        end

        for j = 1:nNodeDoF
            dof = (nodeID-1)*nNodeDoF + j;
            val = bcVals(j);

            if isPrescribed(dof)
                if abs(bcData.prescribedValues(dof) - val) > 1e-12
                    error('Conflicting essential BC detected at global DoF %d.', dof);
                end
            else
                bcData.prescribedValues(dof) = val;
                bcData.prescribedDoFs = [bcData.prescribedDoFs, dof];
                isPrescribed(dof) = true;
            end
        end
    end
end

bcData.prescribedDoFs = unique(bcData.prescribedDoFs);

end