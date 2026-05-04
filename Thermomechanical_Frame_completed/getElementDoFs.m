function eleDoFs = getElementDoFs(eleNodeIDs, nNodeDoF)

nEleNodes = numel(eleNodeIDs);
eleDoFs = zeros(nEleNodes*nNodeDoF,1);

for a = 1:nEleNodes
    for i = 1:nNodeDoF
        eleDoFs((a-1)*nNodeDoF+i) = (eleNodeIDs(a)-1)*nNodeDoF + i;
    end
end

end