function dofInfo = getDofInfo(model)
% Determine nodal DoF layout from model.fieldNames

if ~isfield(model,'fieldNames')
    error('model.fieldNames must be defined.');
end

labels = model.fieldNames(:)';

dofInfo.labels = labels;
dofInfo.nDoFPerNode = numel(labels);

for i = 1:numel(labels)
    dofInfo.map.(labels{i}) = i;
end

end