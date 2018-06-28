% st = struct
% st.faceID=5
% st.cellID_old=G.faces.neighbors(st.faceID,1)
% st.cellID_new=G.faces.neighbors(st.faceID,2)
% st.x=G.faces.centroids(st.faceID,:)
% st.tau=0

function streamline(G, porosity, flux, st)

%-----------------%
% 1. find the barycentric coordinates of x_0 = x(tau_0) of teh entry point wrt
%E.
nod=nodesofcell(G, st.cellID_new);
a=G.nodes.coords(nod(1),:);
b=G.nodes.coords(nod(2),:);
c=G.nodes.coords(nod(3),:);
d=G.nodes.coords(nod(4),:);
lambda = ToBarycentric(a,b,c,d, st.x);

%-----------------%
% 2. compute the velocity in barycentric coordinates
velocity = getVelocity(G, porosity, flux, st.cellID_new);

%-----------------%
% 3.a ffor all i such taht v_i<0, copmture the intersection times
% tau_i = - x_i(tau_0)/v)i
taunew=-lambda./velocity;
% For all other i set tau_i = -Inf
taunew(velocity>=0)=-Inf;

%-----------------%
% 3.b find minimum non-negative intersction time tau_m
validentries=find(taunew>0);
[taunew,ind]=min(taunew(validentries));
indices=validentries(ind);
%todo treat multiple indices

%-----------------%
% 3.c the exit point is given in barycentric coordinates wrt the cell by
% x(tau_m) = x_0 +tau_m*v
lambda_new=lambda+taunew(ind)


end

function lambda = ToBarycentric(a,b,c,d, x)
A=[[a;b;c;d]';[1,1,1,1]];
x=[x';1];
lambda = A\x;
end
function x = FromBarycentric(a,b,c,d, lambda)
A=[[a;b;c;d]';[1,1,1,1]];
x=A*lambda;
x=x(1:3)';
end



function velocity = getVelocity(G, porosity, flux, cellID)
%normal = G.faces.normals(faceID,:);
%neighbors = G.faces.neighbors(faceID,:);
flux=sum(flux(nodesofcell(G, cellID),:),2);
%%if neighbors(1)==0 -> no neighbor
%if neighbors(1) ~= cellID && neighbors(2) ~= cellID
%    return 1;
%end
velocity=-flux/(3*G.cells.volumes(cellID)*porosity(cellID));
%v1=-1/(3*G.cells.volumes(neighbors(1))*porosity(neighbors(1)))*flux
%v2=-1/(3*G.cells.volumes(neighbors(2))*porosity(neighbors(2)))*(-flux)
end