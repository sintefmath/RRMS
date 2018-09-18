function st = streamline(G, porosity, flux, x, faceID, stopcellID, maxTof)
    totflux = sum(flux,2);
    assert(length(totflux)==G.faces.num, "We need a flux per face, aborting!");


    tau=0;

    %%%
    % at each step we have a faceid, cellid
    %%%
    st = struct;
    st.LAMBDA = [];
    st.X = [];
    st.TOF = [];
    st.PARSEDCELLS = [];
    st.PARSEDFACES = [];
    st.breakreason = 'ok';

    while true



        %-----------------%
        % 2. compute the velocity in barycentric coordinates
        neighbors = G.faces.neighbors(faceID,:);
        cellID = neighbors(1+(totflux(faceID)>=0));%if flux is positive, it flows to cell neighbors(2), otherwise neighbors(1).
        %%% stopping criterion, on boundary or marked cell or max tof
        if min(neighbors) == 0
            st.breakreason='boundary cell';
            break;
        end
        if ismember(cellID,stopcellID)
            st.breakreason='well cell';
            break;
        end
        if tau>=maxTof
            st.breakreason='max time reached';
            break;
        end
        facesofcell = get_faces_of_cell(G, cellID);
        neighborsoffaces = G.faces.neighbors(facesofcell,:);
        signu=neighborsoffaces(:,1)==cellID;

        %-----------------%
        % 1. find the barycentric coordinates of x_0 = x(tau_0) of the entry
        % point wrt E.
        nodesofcell = get_nodes_of_cell(G, facesofcell);
        A = G.nodes.coords(nodesofcell(:),:);

        nodesface1=get_nodes_of_face(G, facesofcell(1));
        nodesface2=get_nodes_of_face(G, facesofcell(2));
        nodesface3=get_nodes_of_face(G, facesofcell(3));
        nodesface4=get_nodes_of_face(G, facesofcell(4));
        inda=~ismember(nodesofcell,nodesface1);
        a=A(inda,:);
        indb=~ismember(nodesofcell,nodesface2);
        b=A(indb,:);
        indc=~ismember(nodesofcell,nodesface3);
        c=A(indc,:);
        indd=~ismember(nodesofcell,nodesface4);
        d=A(indd,:);
        lambda = ToBarycentric(a,b,c,d, x);

        %-----------------%
        % 3.a for all i such taht v_i<0, compute the intersection times
        % tau_i = - x_i(tau_0)/v_i
        taunew = -lambda./velocity;
        % For all other i set tau_i = -Inf
        taunew(velocity>=0) = -Inf;
        %taunew(find(lambda<1e-10)) = -Inf;

        %-----------------%
        % 3.b find minimum non-negative intersection time tau_m
        validentries = find(taunew>0);
        [taunew,ind] = min(taunew(validentries));
        indices = validentries(ind);
        %TODO: treat multiple indices
        if length(indices)>1
            st.breakreason='multiple indices found';
            disp(st.breakreason);
            break;
        end
        if length(indices)==0
            st.breakreason='no index found';
            disp(st.breakreason);
            break;
        end

        st.LAMBDA(:,end+1) = lambda;
        st.X(:,end+1) = x;
        st.TOF(end+1) = tau;
        st.PARSEDCELLS(end+1) = cellID;
        st.PARSEDFACES(end+1) = faceID;

        %-----------------%
        % 3.c the exit point is given in barycentric coordinates wrt the cell by
        % x(tau_m) = x_0 +tau_m*v
        tau = tau + taunew;
        faceID = facesofcell(indices);
        lambda = lambda + taunew*velocity;
        x = FromBarycentric(a,b,c,d, lambda);

    end

end

