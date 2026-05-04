function nodalField = averageFieldToNodes(model, elementResults, fieldName)
% Average element result field to nodes
%
% Input:
%   model
%   elementResults : output of postprocessModel
%   fieldName      : name of element result field, e.g.
%                    'T', 'vonMises', 'stress', 'strain',
%                    'thermalStrain', 'heatFlux'
%
% Output:
%   nodalField : numNodes x nComp matrix
%
% Examples:
%   Tnodal     = averageFieldToNodes(model, results, 'T');          % numNodes x 1
%   vmNodal    = averageFieldToNodes(model, results, 'vonMises');   % numNodes x 1
%   stressN    = averageFieldToNodes(model, results, 'stress');     % numNodes x 3
%   strainN    = averageFieldToNodes(model, results, 'strain');     % numNodes x 3
%   heatFluxN  = averageFieldToNodes(model, results, 'heatFlux');   % numNodes x 2

numNodes = size(model.nodesInfo,1);

% Determine number of components from first valid element result
nComp = [];
for e = 1:numel(elementResults)
    if isfield(elementResults(e).value, fieldName)
        val = elementResults(e).value.(fieldName);
        val = val(:);
        nComp = numel(val);
        break;
    end
end

if isempty(nComp)
    error('Field "%s" was not found in elementResults.', fieldName);
end

sumVals = zeros(numNodes, nComp);
countVals = zeros(numNodes, 1);

for e = 1:numel(elementResults)
    if ~isfield(elementResults(e).value, fieldName)
        continue;
    end

    nodeIDs = elementResults(e).nodeIDs;
    val = elementResults(e).value.(fieldName);
    val = val(:)';   % row vector: 1 x nComp

    if numel(val) ~= nComp
        error('Inconsistent component count in field "%s".', fieldName);
    end

    for a = 1:numel(nodeIDs)
        nodeID = nodeIDs(a);
        idx = find(model.nodesInfo(:,1) == nodeID, 1);

        if isempty(idx)
            error('Node ID %d was not found in nodesInfo.', nodeID);
        end

        sumVals(idx,:) = sumVals(idx,:) + val;
        countVals(idx) = countVals(idx) + 1;
    end
end

nodalField = zeros(numNodes, nComp);
for i = 1:numNodes
    if countVals(i) > 0
        nodalField(i,:) = sumVals(i,:) / countVals(i);
    end
end

end