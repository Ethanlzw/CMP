function dofInfo = getDofInfo(model)
% getDofInfo  Build nodal DOF layout from model.fieldNames.
%
%   dofInfo.labels      – cell array of field names, e.g. {'ux','uy'}
%   dofInfo.nDofPerNode – number of DOFs per node
%   dofInfo.map         – struct mapping label -> local DOF index

if ~isfield(model, 'fieldNames') || isempty(model.fieldNames)
    error('model.fieldNames must be a non-empty cell array of field labels.');
end

labels = model.fieldNames(:)';   % ensure row cell array

dofInfo.labels      = labels;
dofInfo.nDofPerNode = numel(labels);
dofInfo.map         = struct();

for i = 1:numel(labels)
    dofInfo.map.(labels{i}) = i;
end

end
