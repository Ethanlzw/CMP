function [K, F] = buildSystem(model)
% buildSystem  Assemble global FE system and add external nodal forces.
%
%   K – global stiffness matrix
%   F – total load vector (internal body forces + external nodal forces)

[K, F] = assembleSystem(model);

if isfield(model, 'bcData') && isfield(model.bcData, 'forceVector')
    F = F + model.bcData.forceVector;
end

end
