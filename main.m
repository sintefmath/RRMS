%%% Copyright 2018 Equinor ASA


%viscosities=[1 10]
%relperms=[2 2]
function [G, porosity, states] = main(viscosities, relperms, vtkfile)


cd /home/franzf/Downloads/checkout-mrst/mrst-bitbucket/mrst-core/
startup
mrstPath('register', 'agmg','/home/franzf/git_repos/agmg/')
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

filename_sim_steps=strcat('sim_2ph_steps', num2str(relperms(1)), '_', num2str(relperms(2)), '_visc_', num2str(viscosities(1)), '_', num2str(viscosities(2)), '_sim.mat');
figure(1); PwellSols_100_1 = sim_2ph_steps(model, W, viscosities, relperms,100,1);
figure(2); cla; hold on; info_100_1 = infoandplot(PwellSols_100_1);

figure(1); PwellSols_100_10 = sim_2ph_steps(model, W, viscosities, relperms,100,10);
figure(2); info_100_10 = infoandplot(PwellSols_100_10);

figure(1); PwellSols_100_100 = sim_2ph_steps(model, W, viscosities, relperms,100,100);
figure(2); info_100_100 = infoandplot(PwellSols_100_100);


figure(1); PwellSols_10_1 = sim_2ph_steps(model, W, viscosities, relperms,10,1);
figure(2); info_10_1 = infoandplot(PwellSols_10_1);

figure(1); PwellSols_10_10 = sim_2ph_steps(model, W, viscosities, relperms,10,10);
figure(2); info_10_10 = infoandplot(PwellSols_10_10);

figure(1); PwellSols_1_1 = sim_2ph_steps(model, W, viscosities, relperms,1,1);
figure(2); info_1_1 = infoandplot(PwellSols_1_1);

legend( ...
    info_100_1, ...
    info_100_10, ...
    info_100_100, ...
    info_10_1, ...
    info_10_10, ...
    info_1_1 ...
    );

filename_sim=strcat('sim_2ph_', num2str(relperms(1)), '_', num2str(relperms(2)), '_visc_', num2str(viscosities(1)), '_', num2str(viscosities(2)), '_sim.mat');

[wellSols, states, report] = sim_2ph(model, W, viscosities, relperms);

if exist(filename_sim, 'file') ~= 2
    [wellSols, states, report] = sim_2ph(model, W, viscosities, relperms);
    save(filename_sim, 'wellSols', 'states', 'report');
else
    load(filename_sim, 'wellSols', 'states', 'report');
end

diagnostics(wellSols, states, report, model, W, 1,0.05,1000,true,200)
%diagnostics(wellSols, states, report, model, W, 20,1,1000,false)
%diagnostics(wellSols, states, report, model, W, 20,0.1,1000,false)
end

function [info] = infoandplot(PwellSols)
value = abs(PwellSols.qWs)./(abs(PwellSols.qOs)+abs(PwellSols.qWs));
plot(PwellSols.t,value);
info=strcat("transport solvs=",num2str(PwellSols.transportSolvs), "pressure solvs=",num2str(PwellSols.pressureSolvs));
end