%%%Copyright 2018 SINTEF AS

% In this example we will solve an incompressible two-phase oil-water
% problem, which consists of an elliptic pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\lambda K\nabla p,$$
%
% where v is the Darcy velocity (total velocity) and lambda is the
% mobility, which depends on the water saturation S.
%
% The saturation equation (conservation of the water phase) is given as:
%
% $$ \phi \frac{\partial S}{\partial t} +
%     \nabla \cdot (f_w(S) v) = q_w$$
%
% where phi is the rock porosity, f is the Buckley-Leverett fractional flow
% function, and q_w is the water source term.
%
% This is an inspired by /modules/incomp/examples/2ph/incompExampleSAIGUP2ph.m

mrstModule add incomp

%% Decide which linear solver to use
linsolve_p = @mldivide;  % Pressure
linsolve_t = @mldivide;  % Transport (implicit)

[G, pressure, pa, velocity, speed, permeability, porosity] = readVTK('z1tracingtt.vtk');

rock.perm = repmat(permeability.* milli*darcy(), [1,3]);
rock.poro = porosity;

%% Set fluid data
% For the two-phase fluid model, we use values:
%
% * densities: [rho_w, rho_o] = [1000 700] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 5] cP.
gravity off
fluid      = initSimpleFluid('mu' , [   1,   5]*centi*poise     , ...
                             'rho', [1000, 700]*kilogram/meter^3, ...
                             'n'  , [   2,   2]);


%% Add wells
%% find all indices that satisfy

maxdist=10;

%% vertical producers (high pressure)
w_x=0;
w_y=0;

productionInx=[];
for i=1:G.cells.num
    dist = (G.cells.centroids(i,1)-w_x)**2;
    dist+= (G.cells.centroids(i,2)-w_y)**2;
    if (sqrt(dist)<maxdist)
        productionInx= [productionInx; i];
    end
end
%WI=ones(length(productionInx),1);
%W = addWell([], G, rock, productionInx, 'Sign', -1, 'Dir', 'z', 'WI', WI);% , 'Radius', 0.1, 'Comp_i', [0, 1]);
W = addWell([], G, rock, productionInx, ...
             'InnerProduct', 'ip_tpf', ...
            'Type', 'bhp', 'Val', 100.0*barsa, ...
            'Radius', 0.1, 'Comp_i', [0, 1]);


%% vertical injectors (low pressure)

w_x=700;
w_y=400;

injectionInx=[];
for i=1:G.cells.num
    dist = (G.cells.centroids(i,1)-w_x)**2;
    dist+= (G.cells.centroids(i,2)-w_y)**2;
    if (sqrt(dist)<maxdist)
        injectionInx= [injectionInx; i];
    end
end

WI=ones(length(injectionInx),1);
%W = addWell(W, G, rock, injectionInx, 'Sign', 1, 'Dir', 'z', 'WI', WI);% , 'Radius', 0.1, 'Comp_i', [1, 0]);
W = addWell(W, G, rock, injectionInx, ...
             'InnerProduct', 'ip_tpf', ...
            'Type', 'bhp' , 'Val', 500.0*barsa, ...
            'Radius', 0.1, 'Dir', 'z', 'Comp_i', [1, 0]);



%% Compute transmissibilities and init reservoir
trans = computeTrans(G, rock, 'Verbose', true);
rSol  = initState(G, W, 0, [0, 1]);

%% Construct pressure and transport solvers
solve_press  = @(x) incompTPFA(x, G, trans, fluid, 'wells', W, ...
                                    'LinSolve', linsolve_p);
solve_transp = @(x, dt) ...
   implicitTransport(x, G, dt, rock, fluid, ...
                     'wells', W, 'LinSolve', linsolve_t);

%% Solve initial pressure
% Solve linear system construced from S and W to obtain solution for flow
% and pressure in the reservoir and the wells.
rSol = solve_press(rSol);



myview = struct('vw',   [-110,18],  ...  % view angle
                'zoom', 1.0,        ...  % zoom
                'asp',  [15 15 2],  ...  % data aspect ratio
                'wh',   50,         ...  % well height above reservoir
                'cb',   'horiz'     ...  % colorbar location
                );

figure(1)
clf
title('pressure')
plotCellData(G, rSol.pressure);
plotWell(G, W, 'height', myview.wh, 'color', 'k');
colorbar
axis tight off, view(myview.vw)
drawnow

figure(2)
clf
title('saturation water')
%plotCellData(G, convertTo(rSol.pressure(1:G.cells.num), barsa));
plotCellData(G, rSol.s(:,1));
plotWell(G, W, 'height', myview.wh, 'color', 'k');
colorbar
axis tight off, view(myview.vw)
drawnow


%% Main loop
% In the main loop, we alternate between solving the transport and the flow
% equations. The transport equation is solved using the standard implicit
% single-point upwind scheme with a simple Newton-Raphson nonlinear solver.
T      = 12*year();
dT     = 1*year();
dTplot = 1*year();

% Start the main loop
t  = 0;  plotNo = 1;
while t < T,
   rSol = solve_transp(rSol, dT);

   % Check for inconsistent saturations
   assert(max(rSol.s(:,1)) < 1+eps && min(rSol.s(:,1)) > -eps);

   % Update solution of pressure equation.
   rSol = solve_press(rSol);

    % Increase time and continue if we do not want to plot saturations
   t = t + dT;
   if ( t < plotNo*dTplot && t <T), continue, end

    figure(1)
    clf
    title('pressure')
    plotCellData(G, rSol.pressure);
    plotWell(G, W, 'height', myview.wh, 'color', 'k');
    colorbar
    axis tight off, view(myview.vw)
    drawnow

    figure(2)
    clf
    title('saturation water')
    %plotCellData(G, convertTo(rSol.pressure(1:G.cells.num), barsa));
    plotCellData(G, rSol.s(:,1));
    plotWell(G, W, 'height', myview.wh, 'color', 'k');
    colorbar
    axis tight off, view(myview.vw)
    drawnow
    %clf, title(strcat('pressure(', num2str(convertTo(t,year)),' years)'))
    %%plotCellData(G, rSol.pressure);
    %%plotCellData(G, rSol.s(:,2));
    %plotCellData(G, rSol.s(:,1), find(rSol.s(:,1) > 0.01))
    %%plotCellData(G, convertTo(rSol.pressure(1:G.cells.num), barsa));
    %plotWell(G, W, 'height', myview.wh, 'color', 'k');
    %colorbar
    %axis tight off, view(myview.vw)
    %drawnow

end
%%

