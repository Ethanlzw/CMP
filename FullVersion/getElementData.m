function elementInfo = getElementData(elementName, model)
% Return metadata and function handle for one element type

dofInfo = getDofInfo(model);

switch elementName
    case 'Tri3_Mech'
        elementInfo.nEleNodes = 3;
        elementInfo.func = @Tri3_Mech;
        elementInfo.requiredFields = {'ux','uy'};

    case 'Quad4_Mech'
        elementInfo.nEleNodes = 4;
        elementInfo.func = @Quad4_Mech;
        elementInfo.requiredFields = {'ux','uy'};

    case 'Tri3_Thermal'
        elementInfo.nEleNodes = 3;
        elementInfo.func = @Tri3_Thermal;
        elementInfo.requiredFields = {'T'};

    case 'Quad4_Thermal'
        elementInfo.nEleNodes = 4;
        elementInfo.func = @Quad4_Thermal;
        elementInfo.requiredFields = {'T'};

    case 'Tri3_ThermoMech'
        elementInfo.nEleNodes = 3;
        elementInfo.func = @Tri3_ThermoMech;
        elementInfo.requiredFields = {'ux','uy','T'};

    case 'Quad4_ThermoMech'
        elementInfo.nEleNodes = 4;
        elementInfo.func = @Quad4_ThermoMech;
        elementInfo.requiredFields = {'ux','uy','T'};

    otherwise
        error('Unsupported element name: %s', elementName);
end

elementInfo.nDoFPerNode = dofInfo.nDoFPerNode;
elementInfo.nEleDoFs = elementInfo.nEleNodes * dofInfo.nDoFPerNode;

if ~isequal(elementInfo.requiredFields, dofInfo.labels)
    error('Element %s is incompatible with model.fieldNames.', elementName);
end

end