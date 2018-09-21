function plotstreamline(G,st, mark)
plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1)
hold on
plot3(st.X(1,:), st.X(2,:), st.X(3,:), 'LineWidth',1)
if mark
    for i=1:length(st.X)
        plot3(st.X(1,i), st.X(2,i), st.X(3,i),'.','MarkerSize',10)
        text(st.X(1,i), st.X(2,i), st.X(3,i),int2str(i))
    end
    view(-45,15)
end
end