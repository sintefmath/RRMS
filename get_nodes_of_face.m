%%%Copyright 2018 SINTEF AS
function nodes = get_nodes_of_face(G, faceID)
    nodes = [G.faces.nodes(3*(faceID-1)+1:3*faceID,:)];
end