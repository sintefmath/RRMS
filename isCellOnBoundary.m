%%%Copyright 2018 SINTEF AS
function val = isCellOnBoundary(G,cellID)

val=false;
faceID=G.cells.faces(4*(cellID-1)+1:4*cellID,:);
for ind=1:length(faceID)
    neighbors=G.faces.neighbors(faceID(ind),:);
    if neighbors(1)==0 || neighbors(2)==0
        val=true;
        return;
    end
end
end