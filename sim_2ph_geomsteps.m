%%%Copyright 2018 SINTEF AS

function [PwellSols] = sim_2ph_steps(model, W, viscosities, relperms, totalsteps)

mrstModule add ad-core ad-blackoil agmg ad-fi ad-props linearsolvers
mrstModule add incomp

po    = 300*barsa;
rSol  = initResSol(model.G, po, [0, 1]);

mrstVerbose false

% * densities: [rho_w, rho_o] = [1000 1000] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 10] cP.
fluid = initSimpleFluid('mu' , viscosities .* centi*poise     , ...
    'rho', [1000, 1000] .* kilogram/meter^3, ...
    'n'  , relperms);

%%
%model.AutoDiffBackend = DiagonalAutoDiffBackend();
%model.AutoDiffBackend = SparseAutoDiffBackend();
model.FacilityModel = UniformFacilityModel(model);
% model.FacilityModel = FacilityModel(model);
model.dsMaxAbs = 0.1;
model.dpMaxRel = 0.1;

%lsolve = setupAMGCL(model);
lsolve = CPRSolverAD('ellipticSolver', AMGCLSolverAD());

nls = NonLinearSolver();
%nls.timeStepSelector = IterationCountTimeStepSelector('firstRampupStep', 5*day, 'targetIterationCount', 5, 'maxTimestep', 30*day);
nls.maxIterations = 12;
nls.LinearSolver = lsolve;
nls.useRelaxation = true;

pv = poreVolume(model.G,model.rock);
T  = sum(pv)/W(2).val;

count = 0;
dT = T/(2^totalsteps);
nextPressuresolve_t = 0;


%% Solve initial pressure in reservoir
gravity off

hT  = computeTrans(model.G, model.rock, 'Verbose', true);

rSol = incompTPFA(rSol, model.G, hT, fluid, 'wells', W);
PwellSols.pressureSolvs = 1;

%% Start the main loop
t  = 0;

wellSols = convertIncompWellSols(W,rSol, fluid);
PwellSols.qWs = wellSols{1}(1).qWs;
PwellSols.qOs = wellSols{1}(1).qOs;
PwellSols.t = t;

PwellSols.transportSolvs = 1;

while t < T,
    %   rSol  = explicitTransport(rSol , model.G, dT, model.rock, fluid, 'wells', W);
    rSol = implicitTransport(rSol, model.G, dT, model.rock, fluid, 'wells', W);
    PwellSols.transportSolvs = PwellSols.transportSolvs + 1;
    disp('transport');
    if (t>=nextPressuresolve_t)
        nextPressuresolve_t = nextPressuresolve_t + T/(2^(totalsteps-count));
        count = count + 1;
        rSol = incompTPFA(rSol, model.G, hT, fluid, 'wells', W);
        PwellSols.pressureSolvs = PwellSols.pressureSolvs + 1;
        disp(' --- pressure');
    end
        
    t = min(t + dT,T);
    
    wellSols = convertIncompWellSols(W,rSol, fluid);
    PwellSols.qWs = [PwellSols.qWs; wellSols{1}(1).qWs];
    PwellSols.qOs = [PwellSols.qOs; wellSols{1}(1).qOs];
    PwellSols.t = [PwellSols.t; t];
    
    if ( t <T), continue, end
    
end
