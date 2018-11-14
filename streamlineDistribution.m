%%%Copyright 2018 SINTEF AS
function states = streamlineDistribution(G, porosity, flux, stopCellIDs, maxTof)

parsedfaceIDs=[];
states = struct;
for i=1:G.faces.num
    disp(i)
    if ~ismember(i,parsedfaceIDs)
        st=streamline(G, porosity, flux, i,stopCellIDs, maxTof);
        parsedfaceIDs=[parsedfaceIDs st.PARSEDFACES];
        tmp = size(st.X);
        if (tmp(2)>2)
            states(end+1).LAMBDA=st.LAMBDA;
            states(end+1).X=st.X;
            states(end+1).TOF=st.TOF;
            states(end+1).PARSEDCELLS=st.PARSEDCELLS;
            states(end+1).PARSEDFACES=st.PARSEDFACES;
        end
    end
end
end