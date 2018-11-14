%%%Copyright 2018 SINTEF AS
function faces = get_faces_of_cell(G, cellID)
    faces = [G.cells.faces(G.cells.facePos(cellID):G.cells.facePos(cellID+1)-1)];
end