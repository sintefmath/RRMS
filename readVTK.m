%%%Copyright 2018 SINTEF AS

function [G, permeability, porosity] = readVTK(filename)
%function [G, pressure, pressure_averaged, velocity, speed, permeability, porosity] = readVTK(filename)

%%% 0) header
file = fopen(filename,'r');
if (file==-1)
    error('Can not open file.');
    return;
end

str = fgets(file);
if ~strcmp(str(3:5), 'vtk')
    error('Not valid vtk file.');    
end

str = fgets(file);
str = strrep(str,sprintf('\n'),'');
if ~strcmp(str(1:17), 'Unstructured Grid')
    error('Not an Unstructured Grid');
end
str = fgets(file);
if ~strcmp(str(1:5), 'ASCII')
    error('Not an ASCII file');
end
str = fgets(file);
if ~strcmp(str(1:25), 'DATASET UNSTRUCTURED_GRID')
    error('Not a DATASET UNSTRUCUTRED_GRID');
end

%%% 1) read points
str=readuntil(file, 'POINTS');
num_points = sscanf(str,'%*s %d');

[points,count] = fscanf(file,'%f %f %f', 3*num_points);
if 3*num_points~=count
    error('Number of points inconsistent.');
end
points = reshape(points, 3, num_points).';

%%% 2) read cells
str=readuntil(file, 'CELLS');
num_cells=sscanf(str,'%*s %d');
num_cell_entries=sscanf(str,'%*s %*d %d');

if 5*num_cells~=num_cell_entries
    error('Number of cells inconsistent.');
end
    
[cells,count] = fscanf(file,'%f %f %f', num_cell_entries);
cells = reshape(cells, 5, num_cells).';

cells = cells(:,2:end)+1;

%%% Create Geometry
G = tetrahedralGrid(points,cells);
G = computeGeometry(G);

% %%% 3) Pressure
% str=readuntil(file, 'SCALARS Corrected_PressureBar');
% str = fgets(file);% reads: LOOKUP_TABLE default
% [pressure,count] = fscanf(file,'%f', num_points);
% %%% 4) Velocity
% str=readuntil(file, 'VECTORS Velocity');
% [velocity,count] = fscanf(file,'%f', [3, num_cell_entries]);
% velocity=velocity.';

%%% 5) Permeability
str=readuntil(file, 'SCALARS Permeability');
str = fgets(file);% reads: LOOKUP_TABLE default
[permeability,count] = fscanf(file,'%f', num_cell_entries);

%%% 5) Porosity
str=readuntil(file, 'SCALARS Porosity');
str = fgets(file);% reads: LOOKUP_TABLE default
[porosity,count] = fscanf(file,'%f', num_cell_entries);

% pressure_averaged=zeros(num_cells,1);
% for i = 1:num_cells
%     ind1=cells(i,1);
%     ind2=cells(i,2);
%     ind3=cells(i,3);
%     ind4=cells(i,4);
%     pressure_averaged(i)=(pressure(ind1)+pressure(ind2)+pressure(ind3)+pressure(ind4))/4.;
% end
% 
% speed=vecnorm(velocity.').';

end

function str = readuntil(file,key)
    str="";
    while true
        if length(str)>=length(key)
            if strcmp(str(1:length(key)), key)
                break
            end
        end
        str = fgets(file);
    end
end

