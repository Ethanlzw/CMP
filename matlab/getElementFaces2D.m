function faces = getElementFaces2D(eleName)
% Return local edge definitions for one 2D element

switch eleName
    case {'Tri3_Mech','Tri3_Thermal','Tri3_ThermoMech','Tri3_Piezo'}
        faces = [1 2;
                 2 3;
                 3 1];

    case {'Quad4_Mech','Quad4_Thermal','Quad4_ThermoMech','Quad4_Piezo'}
        faces = [1 2;
                 2 3;
                 3 4;
                 4 1];

    otherwise
        error('Unsupported element name in getElementFaces2D.');
end

end