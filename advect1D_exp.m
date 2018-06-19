%%% Copyright 2018 Equinor ASA

function [s, wcut, t] = advect1D_exp(s0, x, model, T, varargin)
opt = struct('src', 1, 'tstep' , T/1000, 'tol', 1e-5, 'p0', 300*barsa, 'qp', []);
opt = merge_options(opt, varargin{:});
% 1d advection (e.g., along tof) for flux function f from t=1 to T

[f, df] = getSmoothedFracFlow(model, opt.p0);

n = numel(s0);
assert(n==numel(x)-1);

nstep = round(T/opt.tstep);
tstep = T/nstep;

Df = @(s)spdiags(df(s), 0, n, n);
Dt = spdiags(diff(x(:))/tstep, 0, n, n);
Dx = spdiags(ones(n,1)*[1 -1], [-1 0], n,n);
q  = opt.src*eye(n,1);

r1 = @(s)Dt*s - Dx*f(s) - q;
J = @(s)Dt - Dx*Df(s);

s     = s0;
wcut  = zeros(nstep,1);
t = (1:nstep)*T/nstep;
for k = 1:nstep
    s = s + Dt^-1*(Dx*f(s) + q);
    if ~isempty(opt.qp)
        fs   = f(s);
        wcut(k) = sum(fs.*opt.qp)/sum(opt.qp);
    end
end
end


    


