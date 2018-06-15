%%% Copyright 2018 Equinor ASA

function diagnostics(G, rock, W, wellSols, states, report, model, solnr, binsize, tstep, plotref)

mrstModule add ad-core ad-blackoil deckformat agmg ad-fi ad-props linearsolvers diagnostics mrst-gui

%% DIAGNOSTICS
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

bin_edges   = (0:binsize:100)*year;
bin_centers = .5*(bin_edges(1:end-1)+bin_edges(2:end));
[n,bin_tof] = histc(ttof, bin_edges);
[n,bin_fa] = histc(tfa, bin_edges);

flux_tof = pv./ttof;
flux_fa  = pv./tfa;

binflux_tof = accumarray(bin_tof(bin_tof>0), flux_tof(bin_tof>0))/fac_tof;
binflux_fa = accumarray(bin_fa(bin_fa>0), flux_fa(bin_fa>0))/fac_fa;
% plot tof-distribution
figure(1), hold on
stairs(bin_centers/year, [binflux_tof, binflux_fa]);

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
finalT=20*year;
[s, wcut_tof, t] = advect1D(zeros(nx,1), bin_edges, model, finalT, 'qp', binflux_tof, 'tstep', finalT/tstep);
[s, wcut_fa, t] = advect1D(zeros(nx,1), bin_edges, model, finalT, 'qp', binflux_fa, 'tstep', finalT/tstep);


figure(2), hold on
if (plotref)
    plot(t_sim/year, wcut_sim,'DisplayName','simulation')
end
plot(t/year, wcut_tof, 'DisplayName', strcat('tof, binsize=', num2str(binsize),' solnr=', num2str(solnr)))
plot(t/year, wcut_fa, 'DisplayName', strcat('fa, binsize=', num2str(binsize),' solnr=', num2str(solnr)))

axis([0 20 0 1])
legend;


% [muW, muO] = deal(fluid.muW(po), fluid.muO(po));
% lamW = @(sw)fluid.krW(sw)./muW;
% lamO = @(so)fluid.krO(so)./muO;
% fracflow = @(sw)lamW(sw)./(lamW(sw)+lamO(1-sw));

