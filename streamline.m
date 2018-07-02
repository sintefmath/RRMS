%streamline(G, porosity, states{1}.flux, 124)
function [X, TOF] = streamline(G, porosity, flux, faceID)
totflux = sum(flux,2);
totflux = flux(:,2);
assert(length(totflux)==G.faces.num, "We need a flux per face, aborting!");

%start at the middle of the face
x = G.faces.centroids(faceID,:);
tau=0;

%%%
% at each step we have a faceid, cellid
%%%
st = struct;
st.LAMBDA = [];
st.X = [];
st.TOF = [];
st.PARSEDCELLS = [];

while true
    
    %%% stopping criterion, boundary
    if min(G.faces.neighbors(faceID,:)) == 0
        break;
    end
    
    %-----------------%
    % 2. compute the velocity in barycentric coordinates
    %[velocity, cellID, facesofcell] = nextCell(G, porosity, faceID, totflux);
    neighbors = G.faces.neighbors(faceID,:);
    cellID = neighbors(1+(totflux(faceID)>=0));%if flux is positive, it flows to cell neighbors(2), otherwise neighbors(1).
    facesofcell = [G.cells.faces(G.cells.facePos(cellID):G.cells.facePos(cellID+1)-1)];
    velocity = -totflux(facesofcell)/(3*G.cells.volumes(cellID)*porosity(cellID));
    
    %-----------------%
    % 1. find the barycentric coordinates of x_0 = x(tau_0) of the entry
    % point wrt E.
    nodesofcell = unique([
    [G.faces.nodes(3*(facesofcell(1)-1)+1:3*facesofcell(1),:)];
    [G.faces.nodes(3*(facesofcell(2)-1)+1:3*facesofcell(2),:)];
    [G.faces.nodes(3*(facesofcell(3)-1)+1:3*facesofcell(3),:)];
    [G.faces.nodes(3*(facesofcell(4)-1)+1:3*facesofcell(4),:)];
    ]);
    %a = G.nodes.coords(nodesofcell(1),:);
    %b = G.nodes.coords(nodesofcell(2),:);
    %c = G.nodes.coords(nodesofcell(3),:);
    %d = G.nodes.coords(nodesofcell(4),:);
    A = G.nodes.coords(nodesofcell(:),:);
    
    %TODO: sort such that faces are on opposite sides of missing node
    nodesface1=nodesofface(G, facesofcell(1));
    nodesface2=nodesofface(G, facesofcell(2));
    nodesface3=nodesofface(G, facesofcell(3));
    nodesface4=nodesofface(G, facesofcell(4));
    inda=~ismember(nodesofcell,nodesface1);
    a=A(inda,:);
    indb=~ismember(nodesofcell,nodesface2);
    b=A(indb,:);
    indc=~ismember(nodesofcell,nodesface3);
    c=A(indc,:);
    indd=~ismember(nodesofcell,nodesface4);
    d=A(indd,:);
    lambda = ToBarycentric(a,b,c,d, x);
    
    st.LAMBDA(:,end+1) = lambda;
    planenumind=find(lambda<1e-10);
    st.X(:,end+1) = x;
    st.TOF(end+1) = tau;
    st.PARSEDCELLS(end+1) = cellID;

    %-----------------%
    % 3.a for all i such taht v_i<0, compute the intersection times
    % tau_i = - x_i(tau_0)/v_i
    taunew = -lambda./velocity;
    % For all other i set tau_i = -Inf
    taunew(velocity>=0) = -Inf;
    taunew(planenumind) = -Inf;

    %-----------------%
    % 3.b find minimum non-negative intersection time tau_m
    validentries = find(taunew>0);
    [taunew,ind] = min(taunew(validentries));
    indices = validentries(ind);
    %TODO: treat multiple indices
    if length(indices)>1
        a=1;
        facesofcell(indices)
    end
    if length(indices)==0
        break;
    end

    %-----------------%
    % 3.c the exit point is given in barycentric coordinates wrt the cell by
    % x(tau_m) = x_0 +tau_m*v
    tau = tau + taunew;
    faceID = facesofcell(indices)
    lambda = lambda + taunew*velocity;
    x = FromBarycentric(a,b,c,d, lambda);

end

end

function lambda = ToBarycentric(a,b,c,d, x)
A = [[a;b;c;d]';[1,1,1,1]];
x = [x';1];
lambda = A\x;
end
function x = FromBarycentric(a,b,c,d, lambda)
A = [[a;b;c;d]';[1,1,1,1]];
x = A*lambda;
x = x(1:3)';
end
% function lambda = ToBarycentric(a,b,c, x)
% A = [[a;b;c]';[1,1,1]];
% x = [x';1];
% lambda = A\x;
% end
% function x = FromBarycentric(a,b,c, lambda)
% A = [[a;b;c]';[1,1,1]];
% x = A*lambda;
% x = x(1:2)';
% end

function ret = nodesofface(G, c)
ret = [G.faces.nodes(3*(c-1)+1:3*c,:)];
end