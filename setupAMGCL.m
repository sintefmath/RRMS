function lsolve = setupAMGCL(model)
ncomp = model.water + model.oil + model.gas;

lsolve = AMGCL_CPRSolverAD('block_size', 2, 'maxIterations', 150, 'tolerance', 1e-2);
    lsolve.amgcl_setup.use_drs = true;
    % Optimal?
    lsolve.setCoarsening('aggregation')
%     lsolve.setCoarsening('smoothed_aggregation')
%     lsolve.setCoarsening('ruge_stuben')
%     lsolve.setRelaxation('spai0');
    lsolve.setRelaxation('ilu0');

    lsolve.setSRelaxation('ilu0');
%     lsolve.setSRelaxation('spai0');
    lsolve.amgcl_setup.verbose = false;

%     lsolve.setSRelaxation('ilu0');
    
    lsolve.amgcl_setup.npre = 1;
    lsolve.amgcl_setup.npost = 1;
    lsolve.amgcl_setup.ncycle = 1;
%     lsolve.amgcl_setup.max_levels = 4;
%     lsolve.amgcl_setup.pre_cycles = 2;
    
    lsolve.amgcl_setup.direct_coarse = false;
%     lsolve.amgcl_setup.coarse_enough = 5000;
%     lsolve.setRelaxation('damped_jacobi');

    % Experimental
%    lsolve.setSolver('idrs');
     %lsolve.setSolver('bicgstabl');
     %lsolve.amgcl_setup.bicgstabl_l = 10;
%     lsolve.amgcl_setup.bicgstabl_delta = 0;
%     lsolve.amgcl_setup.bicgstabl_convex = false;

%     lsolve.setSolver('fgmres');
%     lsolve.amgcl_setup.idrs_omega = 0;
%     lsolve.amgcl_setup.idrs_s = 4;
%     lsolve.setRelaxation('ilu0');
%     lsolve.setCoarsening('smoothed_aggregation');
%     lsolve.amgcl_setup.npre = 1;
    lsolve.amgcl_setup.drs_eps_dd = 0.3;
    lsolve.amgcl_setup.drs_eps_ps = 0.1*lsolve.amgcl_setup.drs_eps_dd;
%     lsolve.applyRightDiagonalScaling = true;

    
    
    lsolve.doApplyScalingCPR = true;
    lsolve.trueIMPES = false;
    lsolve.amgcl_setup.use_drs = false;

if isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend')
    %lsolve.reduceToCell = true;
    lsolve.keepNumber = model.G.cells.num*ncomp;
    lsolve.amgcl_setup.active_rows = model.G.cells.num*ncomp;
end

lsolve.verbose = true;
lsolve.amgcl_setup.verbose = false;

end