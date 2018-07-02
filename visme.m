function visme(a,b,c,d,x)
xyz=[a' b' c' d']';
tri=delaunay(xyz);
tetramesh(tri,xyz);
hold on;
plot3(x(1),x(2),x(3),'.','MarkerSize',50);
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
grid on