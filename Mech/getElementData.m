function eleData = getElementData(eleName, model)
% getElementData  Registry: return metadata and function handle for one element type.
%
%   Extend this function to register new element types.

dofInfo = getDofInfo(model);

switch eleName
    case 'Tri3_Mech'
        eleData.nEleNodes      = 3;
        eleData.func           = @Tri3_Mech;
        eleData.requiredFields = {'ux', 'uy'};

    case 'Quad4_Mech'
        eleData.nEleNodes      = 4;
        eleData.func           = @Quad4_Mech;
        eleData.requiredFields = {'ux', 'uy'};

    otherwise
        error('getElementData: unsupported element type ''%s''.', eleName);
end

% Verify element is compatible with model DOF layout
if ~isequal(sort(eleData.requiredFields), sort(dofInfo.labels))
    error('Element ''%s'' requires fields {%s} but model has {%s}.', ...
        eleName, ...
        strjoin(eleData.requiredFields, ', '), ...
        strjoin(dofInfo.labels, ', '));
end

eleData.nDofPerNode = dofInfo.nDofPerNode;
eleData.nEleDofs    = eleData.nEleNodes * dofInfo.nDofPerNode;

end
