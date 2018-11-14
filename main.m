%%%Copyright 2018 SINTEF AS


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

figure(1); cla; hold on;
figure(2); cla; hold on;

if exist('PwellSols_1000_1.mat', 'file') ~= 2
    PwellSols_1000_1 = sim_2ph_steps(model, W, viscosities, relperms,1000,1);
    save('PwellSols_1000_1.mat','PwellSols_1000_1');
else
    load('PwellSols_1000_1.mat','PwellSols_1000_1');
end
info_1000_1 = infoandplot(PwellSols_1000_1, false);

if exist('PwellSols_1000_100.mat', 'file') ~= 2
    PwellSols_1000_100 = sim_2ph_steps(model, W, viscosities, relperms,1000,100);
    save('PwellSols_1000_100.mat','PwellSols_1000_100');
else
    load('PwellSols_1000_100.mat','PwellSols_1000_100');
end
info_1000_100 = infoandplot(PwellSols_1000_100, false);

% PwellSols_100_100 = sim_2ph_steps(model, W, viscosities, relperms,100,100);
% info_100_100 = infoandplot(PwellSols_100_100, false);


% if exist('PwellSols_10_1.mat', 'file') ~= 2
%     PwellSols_10_1 = sim_2ph_steps(model, W, viscosities, relperms,10,1);
%     save('PwellSols_10_1.mat','PwellSols_10_1');
% else
%     load('PwellSols_10_1.mat','PwellSols_10_1');
% end
% info_10_1 = infoandplot(PwellSols_10_1, false);

% PwellSols_10_10 = sim_2ph_steps(model, W, viscosities, relperms,10,10);
% info_10_10 = infoandplot(PwellSols_10_10, false);

% PwellSols_1_1 = sim_2ph_steps(model, W, viscosities, relperms,1,1);
% info_1_1 = infoandplot(PwellSols_1_1, false);

if exist('PwellSols_geom_10.mat', 'file') ~= 2
    PwellSols_geom_10 = sim_2ph_geomsteps(model, W, viscosities, relperms,10);
    save('PwellSols_geom_10.mat','PwellSols_geom_10');
else
    load('PwellSols_geom_10.mat','PwellSols_geom_10');
end
info_geom_10 = infoandplot(PwellSols_geom_10, true);

infostr = strcat('relperm = [',num2str(relperms(1)), ', ', num2str(relperms(2)), ...
                 '], viscosities = [', num2str(viscosities(1)), ', ', num2str(viscosities(2)),']');

figure(1);legend( ...
    info_1000_1, ...
    info_1000_100, ...
    ... %info_100_100, ...
    ... %info_10_1, ...
    ... %info_10_10, ...
    ... %info_1_1, ...
    info_geom_10 ...
    );
xlabel('years');
title(strcat('Water cut, ',infostr));

figure(2);legend( ...
    info_1000_1, ...
    info_1000_100, ...
    ... %info_100_100, ...
    ... %info_10_1, ...
    ... %info_10_10, ...
    ... %info_1_1, ...
    info_geom_10 ...
    );
xlabel('years');
title(strcat('Accumulative oil production, ',infostr));

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

function [info] = infoandplot(PwellSols,geom)
value = abs(PwellSols.qWs)./(abs(PwellSols.qOs)+abs(PwellSols.qWs));
figure(1); plot(PwellSols.t/year(),value);

accumval = (PwellSols.t(2:end)-PwellSols.t(1:end-1)).*abs(PwellSols.qOs(1:end-1));
figure(2); plot(PwellSols.t(2:end)/year(),cumsum(accumval));

info=strcat("#transport solvs=",num2str(PwellSols.transportSolvs), ", #pressure solvs=",num2str(PwellSols.pressureSolvs));
if geom
    info=strcat(info," geom");
end
end
