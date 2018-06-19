%%% Copyright 2018 Equinor ASA

function diagnostics(wellSols, states, report, model, W, solnr, binsize, tstep, plotref)

mrstModule add diagnostics

%% DIAGNOSTICS
D = computeTOFandTracerFirstArrival(states{solnr}, model.G, model.rock, 'wells', W, 'computeWellTOFs', true);
%D = computeTOFandTracer(states{solnr}, model.G, model.rock,  'wells', W);
% get approx tof-distribution
ttof = sum(D.tof,2);
tfa  = D.ifa + D.pfa;
pv = poreVolume(model.G, model.rock);

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

fluxdist = computeDistribution(model, states{solnr}, W, bin_edges);
binflux_pu = fluxdist*(bin_edges(2)-bin_edges(1));

binflux_tof = accumarray(bin_tof(bin_tof>0), flux_tof(bin_tof>0))/fac_tof;
binflux_fa = accumarray(bin_fa(bin_fa>0), flux_fa(bin_fa>0))/fac_fa;
% plot tof-distribution
figure(1)
if (plotref)
    clf;
    hold on;
end
stairs(bin_centers/year, binflux_tof, 'DisplayName', strcat('tof, binsize=', num2str(binsize),' solnr=', num2str(solnr)));
stairs(bin_centers/year, binflux_fa, 'DisplayName', strcat('fa, binsize=', num2str(binsize),' solnr=', num2str(solnr)));
stairs(bin_centers/year, binflux_pu, 'DisplayName', strcat('pu, binsize=', num2str(binsize),' solnr=', num2str(solnr)));
legend;

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
[s_tof_exp, wcut_tof_exp, t_tof_exp] = advect1D_exp(zeros(nx,1), bin_edges, model, finalT, 'qp', binflux_tof, 'tstep', finalT/tstep);
[s_fa_exp, wcut_fa_exp, t_fa_exp] = advect1D_exp(zeros(nx,1), bin_edges, model, finalT, 'qp', binflux_fa, 'tstep', finalT/tstep);
[s_pu_exp, wcut_pu_exp, t_pu_exp] = advect1D_exp(zeros(nx,1), bin_edges, model, finalT, 'qp', binflux_pu, 'tstep', finalT/tstep);
[s_tof_imp, wcut_tof_imp, t_tof_imp] = advect1D_imp(zeros(nx,1), bin_edges, model, finalT, 'qp', binflux_tof, 'tstep', finalT/tstep);
[s_fa_imp, wcut_fa_imp, t_fa_imp] = advect1D_imp(zeros(nx,1), bin_edges, model, finalT, 'qp', binflux_fa, 'tstep', finalT/tstep);
[s_pu_imp, wcut_pu_imp, t_pu_imp] = advect1D_imp(zeros(nx,1), bin_edges, model, finalT, 'qp', binflux_pu, 'tstep', finalT/tstep);


figure(2)
if (plotref)
    clf;
    plot(t_sim/year, wcut_sim,'DisplayName','simulation')
end
hold on;
plot(t_tof_exp/year, wcut_tof_exp, '--', 'DisplayName', strcat('tof exp, binsize=', num2str(binsize),' solnr=', num2str(solnr)))
plot(t_fa_exp/year, wcut_fa_exp, ':', 'DisplayName', strcat('fa exp, binsize=', num2str(binsize),' solnr=', num2str(solnr)))
plot(t_pu_exp/year, wcut_pu_exp, ':', 'DisplayName', strcat('pu exp, binsize=', num2str(binsize),' solnr=', num2str(solnr)))
plot(t_tof_imp/year, wcut_tof_imp, '-.', 'DisplayName', strcat('tof imp, binsize=', num2str(binsize),' solnr=', num2str(solnr)))
plot(t_fa_imp/year, wcut_fa_imp, ':', 'DisplayName', strcat('fa imp, binsize=', num2str(binsize),' solnr=', num2str(solnr)))
plot(t_pu_imp/year, wcut_pu_imp, ':', 'DisplayName', strcat('pu imp, binsize=', num2str(binsize),' solnr=', num2str(solnr)))

axis([0 20 0 1])
legend;


% [muW, muO] = deal(fluid.muW(po), fluid.muO(po));
% lamW = @(sw)fluid.krW(sw)./muW;
% lamO = @(so)fluid.krO(so)./muO;
% fracflow = @(sw)lamW(sw)./(lamW(sw)+lamO(1-sw));

