function face = getElementFaces2D(eleNodeIDs, eleName)
% Return face connectivity for one 2D element

switch eleName
    case {'Tri3_Mech','Tri3_Thermal','Tri3_ThermoMech','Tri3_Piezo'}
        face = eleNodeIDs(1:3);

    case {'Quad4_Mech','Quad4_Thermal','Quad4_ThermoMech','Quad4_Piezo'}
        face = eleNodeIDs(1:4);

    otherwise
        error('Unsupported 2D element type: %s', eleName);
end

end