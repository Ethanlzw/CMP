function plotElementPatches(model, coords, faceColor, edgeColor)
% plotElementPatches  Shared helper: draw all element patches onto current axes.
%
%   coords    – [nNodes × 2] (x,y) coordinates (possibly deformed)
%   faceColor – patch face colour string or RGB
%   edgeColor – patch edge colour string or RGB

hold on;  axis equal;  box on;

for s = 1:size(model.elesSets, 1)
    setField = model.elesSets{s, 1};
    elesSet  = model.elesInfo.(setField);
    for e = 1:size(elesSet, 1)
        en = elesSet(e,:);
        patch(coords(en,1), coords(en,2), faceColor, 'EdgeColor', edgeColor);
    end
end

end
