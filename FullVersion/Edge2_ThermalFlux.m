function [edgeMatrix, edgeVector] = Edge2_ThermalFlux(nodeCoords, edgeData, thickness)
% 2-node boundary element for thermal natural BC
%
% edgeData.type:
%   'flux'      -> prescribed normal heat flux qn
%   'convection'-> h*(T - Tinf)

x1 = nodeCoords(1,1); y1 = nodeCoords(1,2);
x2 = nodeCoords(2,1); y2 = nodeCoords(2,2);
edgeLength = sqrt((x2-x1)^2 + (y2-y1)^2);

edgeMatrix = zeros(2,2);
edgeVector = zeros(2,1);

gaussPoints = [-1/sqrt(3), 1/sqrt(3)];
gaussWeights = [1, 1];

for iGauss = 1:2
    s = gaussPoints(iGauss);
    weight = gaussWeights(iGauss);

    N = [(1-s)/2, (1+s)/2];
    jacobian1D = edgeLength/2;

    if strcmpi(edgeData.type,'flux')
        qn = edgeData.qn;
        edgeVector = edgeVector + N' * qn * jacobian1D * weight * thickness;

    elseif strcmpi(edgeData.type,'convection')
        h = edgeData.h;
        Tinf = edgeData.Tinf;

        edgeMatrix = edgeMatrix + (N' * N) * h * jacobian1D * weight * thickness;
        edgeVector = edgeVector + N' * h * Tinf * jacobian1D * weight * thickness;

    else
        error('Unsupported thermal edgeData.type: %s', edgeData.type);
    end
end

end