%%% Copyright 2018 Equinor ASA

function main()

cd /home/franzf/Downloads/checkout-mrst/mrst-bitbucket/mrst-core/
startup
cd /home/franzf/git_repos/RRMS/
viscosities=[1 10]
relperms=[2 2]

filename=strcat('sim_2ph_relperm_', num2str(relperms(1)), '_', num2str(relperms(2)), '_visc_', num2str(viscosities(1)), '_', num2str(viscosities(2)), '.mat')

if exist(filename, 'file') ~= 2
    [G, pressure, pa, velocity, speed, permeability, porosity] = readVTK('z1tracingtt.vtk');
    sim_2ph(G, permeability, porosity, viscosities, relperms)
end

load(filename, 'wellSols', 'states', 'report', 'rock', 'W', 'model');

diagnostics(G, rock, W, wellSols, states, report, model, 1,1,1000,true)
diagnostics(G, rock, W, wellSols, states, report, model, 20,1,1000,false)
diagnostics(G, rock, W, wellSols, states, report, model, 20,0.1,1000,false)
