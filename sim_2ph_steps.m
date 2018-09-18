%%% Copyright 2018 Equinor ASA

function [wellSols, states, report] = sim_2ph(model, W, viscosities, relperms, transportstepsperpressuresolv)

mrstModule add ad-core ad-blackoil agmg ad-fi ad-props linearsolvers
mrstModule add incomp

%schedule = simpleSchedule(rampupTimesteps(20*year, 60*day, 5), 'W', W);

po    = 300*barsa;
rSol  = initResSol(model.G, po, [0, 1]);
%rSol = incompTPFA(rSol, model.G, hT, fluid, 'wells', W);
mrstVerbose true

% * densities: [rho_w, rho_o] = [1000 700] kg/m^3
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

T      = 300*day();
dT     = T/15;
dTplot = 100*day();  % plot only every 100th day
N      = fix(T/dTplot);
pv     = poreVolume(model.G,model.rock);

%[wellSols, states, report] = simulateScheduleAD(rSol, model, schedule, 'NonLinearSolver', nls);
%rSol = initState(G, W, 0, [0, 1]);

%% Solve initial pressure in reservoir
% Solve linear system construced from S and W to obtain solution for flow
% and pressure in the reservoir and the wells.
gravity off

hT  = computeTrans(model.G, model.rock, 'Verbose', true);

rSol = incompTPFA(rSol, model.G, hT, fluid, 'wells', W);

% Report initial state of reservoir
% subplot(2,1,1), cla
%    plotCellData(model.G, convertTo(rSol.pressure(1:model.G.cells.num), barsa));
%    title('Initial pressure'), view(3)
% 
%subplot(2,1,2), 
cla
   cellNo = rldecode(1:model.G.cells.num, diff(model.G.cells.facePos), 2) .';
   plotCellData(model.G, accumarray(cellNo, ...
      abs(convertTo(faceFlux2cellFlux(model.G, rSol.flux), meter^3/day))));
   title('Initial flux intensity'), view(3)
   drawnow();

%% Transport loop
% We solve the two-phase system using a sequential splitting in which the
% pressure and fluxes are computed by solving the flow equation and then
% held fixed as the saturation is advanced according to the transport
% equation. This procedure is repeated for a given number of time steps
% (here we use 15 equally spaced time steps). The error introduced by this
% splitting of flow and transport can be reduced by iterating each time
% step until e.g., the residual is below a certain user-prescribed
% threshold (this is not done herein).
% T      = 300*day();
% dT     = T/15;
% dTplot = 100*day();  % plot only every 100th day
% N      = fix(T/dTplot);
% pv     = poreVolume(model.G,rock);

%%
% The transport equation will be solved by the single-point upstream method
% with either explicit or implicit time discretizations. Both schemes may
% use internal time steps to obtain a stable discretization. To represent
% the two solutions, we create new solution objects to be used by the
% solver with implicit transport step.
rISol = rSol;

%% Start the main loop
t  = 0; plotNo = 1; hi = 'Implicit: '; he = 'Explicit: ';
e = []; pi = []; pe = [];
pressuresolvs=1;
while t < T,
    if (pressuresolvs >= transportstepsperpressuresolv)
        %   rSol  = explicitTransport(rSol , model.G, dT, model.rock, fluid, 'wells', W);
        rISol = implicitTransport(rISol, model.G, dT, model.rock, fluid, 'wells', W);
        pressuresolvs=0;
    end

   % Check for inconsistent saturations
%   s = [rSol.s(:,1); rISol.s(:,1)];
%   assert(max(s) < 1+eps && min(s) > -eps);

   % Update solution of pressure equation.
%   rSol  = incompTPFA(rSol , model.G, hT, fluid, 'wells', W);
   rISol = incompTPFA(rISol, model.G, hT, fluid, 'wells', W);
   pressuresolvs = pressuresolvs+1;

%    % Measure water saturation in production cells in saturation
%    e = [e; sum(abs(rSol.s(:,1) - rISol.s(:,1)).*pv)/sum(pv)]; %#ok
%    pe = [pe; rSol.s(W(2).cells,1)' ];                 %#ok
%    pi = [pi; rISol.s(W(2).cells,1)'];                 %#ok

   % Increase time and continue if we do not want to plot saturations
   t = t + dT;
   if ( t < plotNo*dTplot && t <T), continue, end

   % Plot saturation
   heading = [num2str(convertTo(t,day)),  ' days'];
   r = 0.01;
%    subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.48]), cla
%    plotCellData(model.G, rSol.s(:,1));
%    view(60,50), axis equal off, title([he heading])
% 
 %subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.48]), 
 cla
   plotCellData(model.G, rISol.s(:,1));
   %view(60,50), 
   axis equal off, title([hi heading])
   drawnow

   plotNo = plotNo+1;
end

% %%
% % As we clearly can see from the plots in the figure, the implicit scheme
% % has much more numerical diffusion than the explicit scheme early in the
% % simulation, but as the time increase, the difference is smaller. To
% % verify this, we can plot the error or the breakthrough curves
% %
% n = size(pe,1);
% pargs = {'MarkerSize',6,'MarkerFaceColor',[.5 .5 .5]};
% subplot(2,1,1),
%    plot(1:n,e*100,'-o', pargs{:}),
%    title('Percentage saturation discrepancy')
% subplot(2,1,2),
%    plot(1:n,pe(:,1),'-o',1:n,pi(:,1),'-s',pargs{:})
%    legend('Explicit','Implicit','Location','NorthWest');
%    title('Water breakthrough at heel'); axis tight


