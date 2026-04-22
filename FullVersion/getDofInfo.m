function dofInfo = getDofInfo(model)
% Determine nodal DoF layout from model.fieldNames
%
% Example:
%   {'ux','uy'}       -> 2 DoFs per node
%   {'T'}             -> 1 DoF per node
%   {'ux','uy','T'}   -> 3 DoFs per node

if ~isfield(model,'fieldNames')
    error('model.fieldNames must be defined.');
end

fieldLabels = model.fieldNames(:)';

nNodeDofs = numel(fieldLabels);

dofInfo = struct();
dofInfo.labels = fieldLabels;
dofInfo.nDoFPerNode = nNodeDofs;

for iField = 1:nNodeDofs
    dofInfo.map.(fieldLabels{iField}) = iField;
end

end