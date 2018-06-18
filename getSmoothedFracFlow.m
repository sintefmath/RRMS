%%% Copyright 2018 Equinor ASA

function varargout = getSmoothedFracFlow(model, p0, ph, ix)
if nargin < 3
    ph = 'W';
end
if nargin < 4
    ix = 1;
end

ap = model.getActivePhases;
switch ph
    case 'W'
        f = model.fluid; 
        p_av  = mean(p0);
        if isfield(model.fluid, 'sWcon')
            sMin = model.fluid.sWcon; 
        else
            sMin = 0;
        end
        if numel(sMin) ~= 1
            sMin = sMin(ix);
        end
        sw  = (sMin:.01:1)';
        
        if ~ap(3) % no gas
            [kw, ko] = model.relPermWO(sw, 1-sw, f);
            [lw, lo] = deal(kw/f.muW(p_av), ko/f.muO(p_av));
        else % assume we are interested in oil/water flow ...
            [kw, ko] = model.relPermWO(sw, 1-sw, f, 'cellInx', repmat(ix, [numel(sw), 1]));
            rs_sat = f.rsSat(p_av, 'cellInx', ix);
            [lw, lo] = deal(kw/f.muW(p_av, 'cellInx', ix), ko/f.muO(p_av, rs_sat, true, 'cellInx', ix));
        end
        f = lw./(lw+lo);
        % fit smooth function of the form a*sc(s)^n/(a*sc(s)^n - (1-sc(s))^n)
        sc   = @(s)max(0, s - sMin)./(1-sMin);
        % coefficent function for given n (minimize misfit by least square)
        coeff = @(n)((f-1).*sc(sw).^n)\(-f.*(1-sc(sw)).^n);
        % smoothed fracflow of samples
        ff    = @(n)(coeff(n)*sc(sw).^n)./(coeff(n)*sc(sw).^n + (1-sc(sw)).^n);
        % find best exponent 1<= n <= 5
        fit   = @(n)norm(ff(n)-f);
        [n, ~, flag] = fminbnd(fit, 1, 5);
        if flag > 0
            a = coeff(n);
            varargout{1} = @(s)(a*sc(s).^n)./(a*sc(s).^n + (max(0,1-sc(s))).^n);
            if nargout > 1
                num  = @(s)a*sc(s).^n;
                dnum = @(s)(a*n/(1-sMin))*sc(s).^(n-1);
                den  = @(s)(a*sc(s).^n + (max(1-sc(s), 0)).^n);
                dden = @(s)(a*n/(1-sMin))*sc(s).^(n-1) - (n/(1-sMin))*(max(0, 1-sc(s))).^(n-1);
                varargout{2} = @(s)(dnum(s).*den(s)-num(s).*dden(s))./(den(s).^2);
            end
        else
            error('Unable to produce frac-flow curves')
            %varargout = getFracFlow(model, state, ph);
        end
    otherwise
        error('Not implemented')
end
end
