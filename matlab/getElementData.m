function eleData = getElementData(eleName, model)
% Return metadata and function handle for one element type

dofInfo = getDofInfo(model);

switch eleName
    case 'Tri3_Mech'
        eleData.nEleNodes = 3;
        eleData.func = @Tri3_Mech;
        eleData.requiredFields = {'ux','uy'};

    case 'Quad4_Mech'
        eleData.nEleNodes = 4;
        eleData.func = @Quad4_Mech;
        eleData.requiredFields = {'ux','uy'};

    case 'Tri3_Thermal'
        eleData.nEleNodes = 3;
        eleData.func = @Tri3_Thermal;
        eleData.requiredFields = {'T'};

    case 'Quad4_Thermal'
        eleData.nEleNodes = 4;
        eleData.func = @Quad4_Thermal;
        eleData.requiredFields = {'T'};

    case 'Tri3_ThermoMech'
        eleData.nEleNodes = 3;
        eleData.func = @Tri3_ThermoMech;
        eleData.requiredFields = {'ux','uy','T'};

    case 'Quad4_ThermoMech'
        eleData.nEleNodes = 4;
        eleData.func = @Quad4_ThermoMech;
        eleData.requiredFields = {'ux','uy','T'};

    case 'Tri3_Piezo'
        eleData.nEleNodes = 3;
        eleData.func = @Tri3_Piezo;
        eleData.requiredFields = {'ux','uy','phi'};

    case 'Quad4_Piezo'
        eleData.nEleNodes = 4;
        eleData.func = @Quad4_Piezo;
        eleData.requiredFields = {'ux','uy','phi'};

    otherwise
        error('Unsupported element name: %s.', eleName);
end

eleData.nDofPerNode = dofInfo.nDofPerNode;
eleData.nEleDofs = eleData.nEleNodes * dofInfo.nDofPerNode;

if ~isequal(eleData.requiredFields, dofInfo.labels)
    error('Element %s is incompatible with model.fieldNames.', eleName);
end

end