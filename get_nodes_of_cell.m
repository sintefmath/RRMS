%%%Copyright 2018 SINTEF AS
function nodes = get_nodes_of_face(G, facesofcell)
    nodes = unique([
        [G.faces.nodes(3*(facesofcell(1)-1)+1:3*facesofcell(1),:)];
        [G.faces.nodes(3*(facesofcell(2)-1)+1:3*facesofcell(2),:)];
        [G.faces.nodes(3*(facesofcell(3)-1)+1:3*facesofcell(3),:)];
        [G.faces.nodes(3*(facesofcell(4)-1)+1:3*facesofcell(4),:)];
        ]);
end