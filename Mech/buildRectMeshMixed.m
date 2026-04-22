function [nodesInfo, elesInfo] = buildRectMeshMixed(L, H, nx, ny, nTriCols)
% buildRectMeshMixed  Structured rectangular mesh with mixed Tri3/Quad4.
%
%   Columns 1..nTriCols  → each cell split into 2 Tri3 triangles.
%   Columns nTriCols+1..nx → Quad4 elements.
%
%   Inputs:
%     L, H     – domain dimensions
%     nx, ny   – number of element columns and rows
%     nTriCols – number of left columns filled with triangles (0 = all Quad4)
%
%   Outputs:
%     nodesInfo – [nNodes × 3]  columns: [nodeId, x, y]
%     elesInfo  – struct with fields setTri [nTri×3] and setQuad [nQuad×4]

dx = L / nx;
dy = H / ny;
nNodesX = nx + 1;
nNodesY = ny + 1;

% ── Node coordinates ──────────────────────────────────────────────────────
%   Nodes are numbered row-by-row (x varies fastest):
%   nodeId = col + row*(nx+1) + 1,  col ∈ [0,nx], row ∈ [0,ny]
nNodes  = nNodesX * nNodesY;
[colIdx, rowIdx] = ndgrid(0:nx, 0:ny);   % colIdx(i,j)=i, rowIdx(i,j)=j
nodeIds   = (1:nNodes)';
nodesInfo = [nodeIds, colIdx(:)*dx, rowIdx(:)*dy];  % [id, x, y]

% ── Connectivity (node index = col + row*(nx+1) + 1) ─────────────────────
nodeIndex = @(col, row) col + row*(nx+1) + 1;   % 0-based col/row

triList  = zeros(2 * nTriCols  * ny, 3);
quadList = zeros((nx - nTriCols) * ny, 4);
iT = 0;  iQ = 0;

for row = 0:ny-1
    for col = 0:nx-1
        n1 = nodeIndex(col,   row  );
        n2 = nodeIndex(col+1, row  );
        n3 = nodeIndex(col+1, row+1);
        n4 = nodeIndex(col,   row+1);

        if col < nTriCols
            triList(iT+1,:) = [n1, n2, n3];
            triList(iT+2,:) = [n1, n3, n4];
            iT = iT + 2;
        else
            quadList(iQ+1,:) = [n1, n2, n3, n4];
            iQ = iQ + 1;
        end
    end
end

elesInfo.setTri  = triList(1:iT, :);
elesInfo.setQuad = quadList(1:iQ, :);

end
