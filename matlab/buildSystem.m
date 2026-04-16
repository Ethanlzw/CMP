function [globalMatrix, globalVector] = buildSystem(model)
% Build FE system including domain and natural BC contributions

[globalMatrix, globalVector] = assembleSystem(model);

if ~isfield(model, 'fieldNames')
    error('model.fieldNames must be defined.');
end

fieldNames = model.fieldNames;

if isequal(fieldNames, {'T'})
    [bcMatrix, bcVector] = assembleNaturalBC_Thermal(model);
    globalMatrix = globalMatrix + bcMatrix;
    globalVector = globalVector + bcVector;

elseif isequal(fieldNames, {'ux','uy'})
    bcVector = assembleNaturalBC_Mech(model);
    globalVector = globalVector + bcVector;

elseif isequal(fieldNames, {'ux','uy','T'})
    [bcMatrix, bcVector] = assembleNaturalBC_ThermoMech(model);
    globalMatrix = globalMatrix + bcMatrix;
    globalVector = globalVector + bcVector;

elseif isequal(fieldNames, {'ux','uy','phi'})
    [bcMatrix, bcVector] = assembleNaturalBC_Piezo(model);
    globalMatrix = globalMatrix + bcMatrix;
    globalVector = globalVector + bcVector;

else
    error('Unsupported model.fieldNames definition.');
end

end