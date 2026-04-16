function nodalValue = averageFieldToNodes(model, eleValues)
% Average element scalar values to nodes

nodesNum = size(model.nodesInfo, 1);
nodalValue = zeros(nodesNum, 1);
counter = zeros(nodesNum, 1);

for setI = 1:numel(model.elesSetName)
    setName = model.elesSetName{setI};
    elesSet = model.elesInfo.(setName);

    if isempty(elesSet) || ~isfield(eleValues, setName)
        continue;
    end

    vals = eleValues.(setName);

    for eleI = 1:size(elesSet,1)
        eleNodes = elesSet(eleI,:);
        v = vals(eleI);

        nodalValue(eleNodes) = nodalValue(eleNodes) + v;
        counter(eleNodes) = counter(eleNodes) + 1;
    end
end

counter(counter == 0) = 1;
nodalValue = nodalValue ./ counter;

end