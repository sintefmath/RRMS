function c=nodesofcell(G, cellID)
%facesofcell=@(c) [G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1)];
%nodesofface=@(c) [G.faces.nodes(3*c-1:3*(c+1)-2,:)];
a=[G.cells.faces(G.cells.facePos(cellID):G.cells.facePos(cellID+1)-1)];%facesofcell(c);
c=unique([
    [G.faces.nodes(3*(a(1)-1)+1:3*a(1),:)];
    [G.faces.nodes(3*(a(2)-1)+1:3*a(2),:)];
    [G.faces.nodes(3*(a(3)-1)+1:3*a(3),:)];
    [G.faces.nodes(3*(a(4)-1)+1:3*a(4),:)];
    ]);
%c=unique([nodesofface(a(1)); nodesofface(a(2)); nodesofface(a(3));nodesofface(a(4));]);
end