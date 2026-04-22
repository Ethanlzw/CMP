function [edgeMatrix, edgeVector] = Edge2_ThermoMech(nodeCoords, edgeData, thickness)
% Boundary element for thermo-mechanical natural BC
% DoF order per node: [ux, uy, T]

x1 = nodeCoords(1,1); y1 = nodeCoords(1,2);
x2 = nodeCoords(2,1); y2 = nodeCoords(2,2);
edgeLength = sqrt((x2-x1)^2 + (y2-y1)^2);

edgeMatrix = zeros(6,6);
edgeVector = zeros(6,1);

gaussPoints = [-1/sqrt(3), 1/sqrt(3)];
gaussWeights = [1, 1];

for iGauss = 1:2
    s = gaussPoints(iGauss);
    weight = gaussWeights(iGauss);

    N = [(1-s)/2, (1+s)/2];
    jacobian1D = edgeLength/2;

    if isfield(edgeData,'traction')
        traction = edgeData.traction(:);
        NmatMech = [N(1) 0    0    N(2) 0    0;
                    0    N(1) 0    0    N(2) 0];
        edgeVector = edgeVector + NmatMech' * traction * jacobian1D * weight * thickness;
    end

    if isfield(edgeData,'thermalType')
        NmatTherm = [N(1) N(2)];
        thermalDofs = [3 6];

        if strcmpi(edgeData.thermalType,'flux')
            qn = edgeData.qn;
            edgeVector(thermalDofs) = edgeVector(thermalDofs) + NmatTherm' * qn * jacobian1D * weight * thickness;

        elseif strcmpi(edgeData.thermalType,'convection')
            h = edgeData.h;
            Tinf = edgeData.Tinf;
            edgeMatrix(thermalDofs,thermalDofs) = edgeMatrix(thermalDofs,thermalDofs) + ...
                (NmatTherm' * NmatTherm) * h * jacobian1D * weight * thickness;
            edgeVector(thermalDofs) = edgeVector(thermalDofs) + ...
                NmatTherm' * h * Tinf * jacobian1D * weight * thickness;
        else
            error('Unsupported thermalType in Edge2_ThermoMech.');
        end
    end
end

end