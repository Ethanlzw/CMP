function results = postprocessModel(model, solution)
% Element-level postprocessing
%
% Supports:
%   {'ux','uy'}       -> stress/strain/vonMises
%   {'T'}             -> element average temperature, thermal gradient
%   {'ux','uy','T'}   -> stress/strain/vonMises + temperature

fieldNames = model.fieldNames;
dofInfo = getDofInfo(model);

nodesInfo = model.nodesInfo;
elementTypeTable = model.elesType;
elementSetNames = model.elesSetName;
nElementSets = size(elementTypeTable,1);

nElementsTotal = 0;
for iSet = 1:nElementSets
    nElementsTotal = nElementsTotal + size(model.elesInfo.(elementSetNames{iSet}),1);
end

results = struct();
results.eleID = zeros(nElementsTotal,1);
results.eleName = cell(nElementsTotal,1);
results.centroid = zeros(nElementsTotal,2);

if isequal(fieldNames, {'ux','uy'}) || isequal(fieldNames, {'ux','uy','T'})
    results.strain = zeros(nElementsTotal,3);
    results.stress = zeros(nElementsTotal,3);
    results.vonMises = zeros(nElementsTotal,1);
end

if isequal(fieldNames, {'T'}) || isequal(fieldNames, {'ux','uy','T'})
    results.temperature = zeros(nElementsTotal,1);
    results.tempGrad = zeros(nElementsTotal,2);
end

if isfield(model,'problem') && isfield(model.problem,'mechType')
    mechType = model.problem.mechType;
else
    mechType = 'planeStress';
end

if isfield(model,'problem') && isfield(model.problem,'deltaTRef')
    deltaTRef = model.problem.deltaTRef;
else
    deltaTRef = 0.0;
end

iResult = 0;
for iSet = 1:nElementSets
    elementName = elementTypeTable{iSet,2};
    elementInfo = getElementData(elementName, model);
    nElementNodes = elementInfo.nEleNodes;
    elementSet = model.elesInfo.(elementSetNames{iSet});

    for iElement = 1:size(elementSet,1)
        iResult = iResult + 1;

        elementID = elementSet(iElement,1);
        elementNodeIDs = elementSet(iElement,2:1+nElementNodes);

        elementNodesInfo = zeros(nElementNodes, size(nodesInfo,2));
        for a = 1:nElementNodes
            nodeMask = nodesInfo(:,1) == elementNodeIDs(a);
            elementNodesInfo(a,:) = nodesInfo(nodeMask,:);
        end

        elementDofs = getElementDof(elementNodeIDs, dofInfo);
        uElement = solution(elementDofs);

        materialMask = model.section(:,1) == elementID;
        material = model.section(materialMask,2:end);

        results.eleID(iResult) = elementID;
        results.eleName{iResult} = elementName;
        results.centroid(iResult,:) = mean(elementNodesInfo(:,2:3),1);

        if isequal(fieldNames, {'ux','uy'})
            [strain, stress, vm] = evaluateMechanicalFields(elementName, elementNodesInfo, uElement, material, mechType);
            results.strain(iResult,:) = strain(:).';
            results.stress(iResult,:) = stress(:).';
            results.vonMises(iResult) = vm;

        elseif isequal(fieldNames, {'T'})
            [temperature, gradT] = evaluateThermalFields(elementName, elementNodesInfo, uElement);
            results.temperature(iResult) = temperature;
            results.tempGrad(iResult,:) = gradT(:).';

        elseif isequal(fieldNames, {'ux','uy','T'})
            [strain, stress, vm, temperature, gradT] = ...
                evaluateThermoMechanicalFields(elementName, elementNodesInfo, uElement, material, mechType, deltaTRef);

            results.strain(iResult,:) = strain(:).';
            results.stress(iResult,:) = stress(:).';
            results.vonMises(iResult) = vm;
            results.temperature(iResult) = temperature;
            results.tempGrad(iResult,:) = gradT(:).';

        else
            error('Unsupported fieldNames in postprocessModel.');
        end
    end
end

end

function [strain, stress, vm] = evaluateMechanicalFields(elementName, elementNodesInfo, uElement, material, mechType)

if strcmp(elementName,'Tri3_Mech')
    B = getTri3Bmech(elementNodesInfo);
    displacementDofs = [1 2 3 4 5 6];
    strain = B * uElement(displacementDofs);

elseif strcmp(elementName,'Quad4_Mech')
    B = getQuad4BmechAtCenter(elementNodesInfo);
    displacementDofs = 1:8;
    strain = B * uElement(displacementDofs);

else
    error('Unsupported elementName in evaluateMechanicalFields.');
end

E = material(1);
nu = material(2);
D = getElasticMatrix2D(E, nu, mechType);
stress = D * strain;
vm = computeVonMises2D(stress, material, mechType);

end

function [temperature, gradT] = evaluateThermalFields(elementName, elementNodesInfo, uElement)

if strcmp(elementName,'Tri3_Thermal')
    Bth = getTri3Bth(elementNodesInfo);
    temperatureDofs = [1 2 3];
    gradT = Bth * uElement(temperatureDofs);
    temperature = mean(uElement(temperatureDofs));

elseif strcmp(elementName,'Quad4_Thermal')
    Bth = getQuad4BthAtCenter(elementNodesInfo);
    temperatureDofs = [1 2 3 4];
    gradT = Bth * uElement(temperatureDofs);
    temperature = mean(uElement(temperatureDofs));

else
    error('Unsupported elementName in evaluateThermalFields.');
end

end

function [strain, stress, vm, temperature, gradT] = evaluateThermoMechanicalFields(elementName, elementNodesInfo, uElement, material, mechType, deltaTRef)

if strcmp(elementName,'Tri3_ThermoMech')
    Bmech = getTri3Bmech(elementNodesInfo);
    Bth = getTri3Bth(elementNodesInfo);
    displacementDofs = [1 2 4 5 7 8];
    temperatureDofs = [3 6 9];

elseif strcmp(elementName,'Quad4_ThermoMech')
    Bmech = getQuad4BmechAtCenter(elementNodesInfo);
    Bth = getQuad4BthAtCenter(elementNodesInfo);
    displacementDofs = [1 2 4 5 7 8 10 11];
    temperatureDofs = [3 6 9 12];

else
    error('Unsupported elementName in evaluateThermoMechanicalFields.');
end

mechanicalDisplacements = uElement(displacementDofs);
temperatureValues = uElement(temperatureDofs);

temperature = mean(temperatureValues);
gradT = Bth * temperatureValues;
strain = Bmech * mechanicalDisplacements;

E = material(1);
nu = material(2);
alpha = material(4);
D = getElasticMatrix2D(E, nu, mechType);
thermalStrain = alpha * (temperature - deltaTRef) * [1; 1; 0];
stress = D * (strain - thermalStrain);
vm = computeVonMises2D(stress, material, mechType);

end

function B = getTri3Bmech(elementNodesInfo)

x1 = elementNodesInfo(1,2); y1 = elementNodesInfo(1,3);
x2 = elementNodesInfo(2,2); y2 = elementNodesInfo(2,3);
x3 = elementNodesInfo(3,2); y3 = elementNodesInfo(3,3);

A = [1 x1 y1;
     1 x2 y2;
     1 x3 y3];
area = det(A) / 2;

b1 = y2 - y3;  c1 = x3 - x2;
b2 = y3 - y1;  c2 = x1 - x3;
b3 = y1 - y2;  c3 = x2 - x1;

B = 1/(2*area) * ...
    [b1  0   b2  0   b3  0;
     0   c1  0   c2  0   c3;
     c1  b1  c2  b2  c3  b3];

end

function Bth = getTri3Bth(elementNodesInfo)

x1 = elementNodesInfo(1,2); y1 = elementNodesInfo(1,3);
x2 = elementNodesInfo(2,2); y2 = elementNodesInfo(2,3);
x3 = elementNodesInfo(3,2); y3 = elementNodesInfo(3,3);

A = [1 x1 y1;
     1 x2 y2;
     1 x3 y3];
area = det(A) / 2;

b1 = y2 - y3;  c1 = x3 - x2;
b2 = y3 - y1;  c2 = x1 - x3;
b3 = y1 - y2;  c3 = x2 - x1;

Bth = 1/(2*area) * [b1 b2 b3;
                    c1 c2 c3];

end

function B = getQuad4BmechAtCenter(elementNodesInfo)

xCoords = elementNodesInfo(:,2);
yCoords = elementNodesInfo(:,3);
s = 0; t = 0;

dNds = 0.25 * [-(1-t),  (1-t), (1+t), -(1+t)];
dNdt = 0.25 * [-(1-s), -(1+s), (1+s),  (1-s)];

J = [dNds; dNdt] * [xCoords yCoords];
dNdxy = J \ [dNds; dNdt];
dNdx = dNdxy(1,:);
dNdy = dNdxy(2,:);

B = [dNdx(1) 0        dNdx(2) 0        dNdx(3) 0        dNdx(4) 0;
     0        dNdy(1) 0        dNdy(2) 0        dNdy(3) 0        dNdy(4);
     dNdy(1) dNdx(1) dNdy(2) dNdx(2) dNdy(3) dNdx(3) dNdy(4) dNdx(4)];

end

function Bth = getQuad4BthAtCenter(elementNodesInfo)

xCoords = elementNodesInfo(:,2);
yCoords = elementNodesInfo(:,3);
s = 0; t = 0;

dNds = 0.25 * [-(1-t),  (1-t), (1+t), -(1+t)];
dNdt = 0.25 * [-(1-s), -(1+s), (1+s),  (1-s)];

J = [dNds; dNdt] * [xCoords yCoords];
Bth = J \ [dNds; dNdt];

end

function vm = computeVonMises2D(stress, material, mechType)

sx = stress(1);
sy = stress(2);
txy = stress(3);

switch lower(mechType)
    case lower('planeStress')
        sz = 0.0;
    case lower('planeStrain')
        nu = material(2);
        sz = nu * (sx + sy);
    otherwise
        error('Unsupported mechType: %s', mechType);
end

vm = sqrt(0.5 * ((sx - sy)^2 + (sy - sz)^2 + (sz - sx)^2) + 3 * txy^2);

end