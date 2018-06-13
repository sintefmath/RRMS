

mrstModule add ad-core ad-blackoil deckformat agmg ad-fi ad-props linearsolvers diagnostics mrst-gui

addpath /Users/steink/Desktop/Documents/data/RRMS/
%%
fluid = initSimpleADIFluid('phases', 'WO', 'n', [2 2], 'mu', [1 1]*centi*poise);

[G, pressure, pa, velocity, speed, permeability, porosity] = readVTK('z1tracingtt.vtk');

rock.perm = repmat(permeability.* milli*darcy(), [1,3]);
rock.poro = porosity;


model = TwoPhaseOilWaterModel(G, rock, fluid);

%% find all indices that satisfy

maxdist=10;

% vertical producers (high pressure)
w_x=0;
w_y=0;

productionInx=[];
for i=1:G.cells.num
    dist = (G.cells.centroids(i,1)-w_x).^2;
    dist = dist + (G.cells.centroids(i,2)-w_y).^2;
    if (sqrt(dist)<maxdist)
        productionInx= [productionInx; i];
    end
end
%WI=ones(length(productionInx),1);
%W = addWell([], G, rock, productionInx, 'Sign', -1, 'Dir', 'z', 'WI', WI);% , 'Radius', 0.1, 'Comp_i', [0, 1]);
W = addWell([], G, rock, productionInx, ...
             'InnerProduct', 'ip_tpf', ...
            'Type', 'bhp', 'Val', 100.0*barsa, ...
            'Radius', 0.1, 'Comp_i', [0, 1], 'Sign', -1, 'Name', 'P');


% vertical injectors (low pressure)

w_x=700;
w_y=400;

injectionInx=[];
for i=1:G.cells.num
    dist = (G.cells.centroids(i,1)-w_x).^2;
    dist = dist + (G.cells.centroids(i,2)-w_y).^2;
    if (sqrt(dist)<maxdist)
        injectionInx= [injectionInx; i];
    end
end

WI=ones(length(injectionInx),1);
%W = addWell(W, G, rock, injectionInx, 'Sign', 1, 'Dir', 'z', 'WI', WI);% , 'Radius', 0.1, 'Comp_i', [1, 0]);
W = addWell(W, G, rock, injectionInx, ...
             'InnerProduct', 'ip_tpf', ...
            'Type', 'rate' , 'Val', 2000*meter^3/day, ...
            'Radius', 0.1, 'Dir', 'z', 'Comp_i', [1, 0], 'Sign', 1, 'Name', 'I');


%%        

schedule = simpleSchedule(rampupTimesteps(20*year, 60*day, 5), 'W', W);

po    = 300*barsa;
rSol  = initResSol(G, po, [0, 1]);
mrstVerbose true

%%
%model.AutoDiffBackend = DiagonalAutoDiffBackend();
%model.AutoDiffBackend = SparseAutoDiffBackend();
model.FacilityModel = UniformFacilityModel(model);
% model.FacilityModel = FacilityModel(model);
model.dsMaxAbs = 0.1;
model.dpMaxRel = 0.1;

%lsolve = setupAMGCL(model);  
lsolve = CPRSolverAD('ellipticSolver', AGMGSolverAD());

nls = NonLinearSolver();
%nls.timeStepSelector = IterationCountTimeStepSelector('firstRampupStep', 5*day, 'targetIterationCount', 5, 'maxTimestep', 30*day);
nls.maxIterations = 12;
nls.LinearSolver = lsolve;
nls.useRelaxation = true;


[wellSols, states, report] = simulateScheduleAD(rSol, model, schedule, 'NonLinearSolver', nls);

%% DIAGNOSTICS
solnr = 20;
D = computeTOFandTracerFirstArrival(states{solnr}, G, rock, 'wells', W, 'computeWellTOFs', true);
%D = computeTOFandTracer(states{solnr}, G, rock,  'wells', W);
% get approx tof-distribution
ttof = sum(D.tof,2);
tfa  = D.ifa + D.pfa;
pv = poreVolume(G, rock);

% adjust ttof to match actual flux
fac_tof = sum(pv./ttof)./wellSols{solnr}(2).qWs;
%ttof = fac_tof*ttof;
fac_fa  = sum(pv./tfa)./wellSols{solnr}(2).qWs;
%tfa = fac_fa*tfa;

bin_edges   = (0:.5:100)*year;
bin_centers = .5*(bin_edges(1:end-1)+bin_edges(2:end));
[n,bin_tof] = histc(ttof, bin_edges);
[n,bin_fa] = histc(tfa, bin_edges);

flux_tof = pv./ttof;
flux_fa  = pv./tfa;

binflux_tof = accumarray(bin_tof(bin_tof>0), flux_tof(bin_tof>0))/fac_tof;
binflux_fa = accumarray(bin_fa(bin_fa>0), flux_fa(bin_fa>0))/fac_fa;
% plot tof-distribution
figure, stairs(bin_centers/year, [binflux_tof, binflux_fa]);

%%
% minimum swat = .2 and maximum swat = .8, so effective porevolume is .6*pv
%wcut_tof = cumsum(binflux_tof)/wellSols{solnr}(2).qWs;
%wcut_fa  = cumsum(binflux_fa)/wellSols{solnr}(2).qWs;
%t_approx    = bin_centers/year;

%wcut_sim = cellfun(@(x)x(1).wcut, wellSols);
qo = cellfun(@(x)x(1).qOs, wellSols);
qw = cellfun(@(x)x(1).qWs, wellSols);
wcut_sim = qw./(qo+qw);
t_sim    = report.ReservoirTime;


nx = numel(bin_centers);
[s, wcut_tof, t] = advect1D(zeros(nx,1), bin_edges, model, 20*year, 'qp', binflux_tof);
[s, wcut_fa, t] = advect1D(zeros(nx,1), bin_edges, model, 20*year, 'qp', binflux_fa);


figure, hold on
plot(t_sim/year, wcut_sim)
plot(t/year, [wcut_tof, wcut_fa])
axis([0 20 0 1])
legend({'simulation', 'diagnostics-tof', 'diagnostics-fa'})


% [muW, muO] = deal(fluid.muW(po), fluid.muO(po));
% lamW = @(sw)fluid.krW(sw)./muW;
% lamO = @(so)fluid.krO(so)./muO;
% fracflow = @(sw)lamW(sw)./(lamW(sw)+lamO(1-sw));

