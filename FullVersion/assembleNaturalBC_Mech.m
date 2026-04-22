function bcVector = assembleNaturalBC_Mech(model)
% Assemble mechanical natural BC contributions from model.natBC
%
% model.natBC rows:
%   { [node1,node2,...], [tx,ty] }

nNodeDofs = model.elesType{1,3};
if nNodeDofs ~= 2
    error('assembleNaturalBC_Mech is only for mechanical models with 2 DoFs per node.');
end

nNodes = size(model.nodesInfo,1);
bcVector = zeros(nNodes*nNodeDofs,1);

if ~isfield(model,'natBC') || isempty(model.natBC)
    return;
end

if isfield(model,'problem') && isfield(model.problem,'thickness')
    thickness = model.problem.thickness;
else
    thickness = 1.0;
end

for iBC = 1:size(model.natBC,1)
    bcNodeIDs = model.natBC{iBC,1};
    tractionValues = model.natBC{iBC,2};

    if numel(tractionValues) ~= 2
        error('Mechanical natural BC values must be [tx, ty].');
    end

    for iEdge = 1:(numel(bcNodeIDs)-1)
        nodeA = bcNodeIDs(iEdge);
        nodeB = bcNodeIDs(iEdge+1);

        edgeInfo = findEdgeOwner(model, nodeA, nodeB);
        if ~edgeInfo.found
            error('Boundary edge (%d,%d) was not found in any element.', nodeA, nodeB);
        end

        idxA = find(model.nodesInfo(:,1) == nodeA, 1);
        idxB = find(model.nodesInfo(:,1) == nodeB, 1);

        nodeCoords = [model.nodesInfo(idxA,2:3);
                      model.nodesInfo(idxB,2:3)];

        edgeVector = Edge2_MechTraction(nodeCoords, tractionValues, thickness);

        edgeDofs = [(nodeA-1)*nNodeDofs + 1;
                    (nodeA-1)*nNodeDofs + 2;
                    (nodeB-1)*nNodeDofs + 1;
                    (nodeB-1)*nNodeDofs + 2];

        bcVector(edgeDofs) = bcVector(edgeDofs) + edgeVector;
    end
end

end