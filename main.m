%%% Copyright 2018 Equinor ASA

%viscosities=[1 10]
%relperms=[2 2]
function main(viscosities, relperms, vtkfile)

cd /home/franzf/Downloads/checkout-mrst/mrst-bitbucket/mrst-core/
startup
cd /home/franzf/git_repos/RRMS/


filename_model=strcat(vtkfile, '_model.mat');
if exist(filename_model, 'file') ~= 2
    [G, permeability, porosity] = readVTK(vtkfile);
    save(filename_model, 'G', 'permeability', 'porosity');
else
    load(filename_model, 'G', 'permeability', 'porosity');
end

mrstModule add ad-props % initSimpleADIFluid
fluid = initSimpleADIFluid('phases', 'WO', 'n', relperms, 'mu', viscosities*centi*poise);
rock.perm = repmat(permeability.* milli*darcy(), [1,3]);
rock.poro = porosity;
mrstModule add ad-core ad-blackoil
model = TwoPhaseOilWaterModel(G, rock, fluid);

W = setupWells(model);

filename_sim=strcat('sim_2ph_', num2str(relperms(1)), '_', num2str(relperms(2)), '_visc_', num2str(viscosities(1)), '_', num2str(viscosities(2)), '_sim.mat');
if exist(filename_sim, 'file') ~= 2
    [wellSols, states, report] = sim_2ph(model, W, viscosities, relperms);
    save(filename_sim, 'wellSols', 'states', 'report');
else
    load(filename_sim, 'wellSols', 'states', 'report');
end

diagnostics(wellSols, states, report, model, W, 1,1,100,true)
%diagnostics(wellSols, states, report, model, W, 20,1,1000,false)
%diagnostics(wellSols, states, report, model, W, 20,0.1,1000,false)
