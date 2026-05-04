function nodal = extractNodalFields(model, solution)
% Extract nodal fields from global solution vector using model.fieldNames

nNodeDoF = model.elesType{1,3};
numNodes = size(model.nodesInfo,1);

if ~isfield(model,'fieldNames')
    error(['model.fieldNames must be defined, e.g. ', ...
           '{''ux'',''uy''}, {''T''}, {''ux'',''uy'',''T''}, or {''ux'',''uy'',''phi''}.']);
end

fieldNames = model.fieldNames;

if numel(fieldNames) ~= nNodeDoF
    error('Number of model.fieldNames must equal nNodeDoF.');
end

nodal = struct();
nodal.nodeIDs = model.nodesInfo(:,1);
nodal.coords = model.nodesInfo(:,2:end);

for i = 1:nNodeDoF
    nodal.(fieldNames{i}) = zeros(numNodes,1);
end

for n = 1:numNodes
    base = (n-1)*nNodeDoF;
    for i = 1:nNodeDoF
        nodal.(fieldNames{i})(n) = solution(base+i);
    end
end

end