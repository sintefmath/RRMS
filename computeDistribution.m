function fluxdist = computeDistribution(model, state, W, bin_edges)
db = diff(bin_edges(:));
dt = db(1);
assert(all( abs(dt-db)/dt < sqrt(eps) ), 'Function assumes equal bin-sizes')

ws = state.wellSol;

q  = sum(vertcat(ws.flux), 2);
wc = vertcat(W.cells);

qInj  = q.*(q>0);
qProd = q.*(q<0);

s  = zeros(model.G.cells.num, 1);
v  = sum(state.flux(model.operators.internalConn, :), 2);
pv = model.operators.pv;

N = model.operators.N;

upIx = zeros(size(v));
upIx(v>=0) = N(v>=0, 1);
upIx(v<0) = N(v<0, 2);

[nf, nc] = deal(numel(v), numel(s));
C  = sparse( [(1:nf)'; (1:nf)'], N, ones(nf,1)*[1 -1], nf, nc);

[qi, qp] = deal(sparse(wc, 1, qInj, nc, 1) ,sparse(wc, 1, qProd, nc, 1));

s0 = full(qi./pv);
s = initVariablesADI(s);

eq = s + (dt./pv).*(C'*(v.*s(upIx))-s.*qp); %-s0
A = eq.jac{1};
A = decomposition(A);

nSteps = numel(db);
fluxdist = zeros(nSteps, 1);
s = s0;
%fprintf('%4.1f %%', 0)
for k = 1:nSteps
    s = A\s;
    fluxdist(k) = -qp'*s;
    %fprintf('\b\b\b\b\b%4.1f%%', 100*k/nSteps)
end
%fprintf('\n')
end

% 
%     
% 
% eq = (C'*(v.*s(upIx).*tr(upIx))-s.*sum(qp,2))./pvr;
% %eq = (C'*(v.*s(upIx))-s.*sum(qp,2))./pvr;
% %eq = (C'*(v.*s(upIx))-s.*sum(qp,2))./pvr;
% s = double(s);
% A  = - eq.jac{1};
% M = (speye(numel(s)) -dt*A);
% 
% injIx  = wsel.injectorIx;
% prodIx = wsel.producerIx;
% if isempty(prodIx) || isempty(injIx)
%     warning('Empty well selection, no distribution ...')
%     [t, qr, pvr, q, pv] = deal([]);
%     return;
% end
% state  = d.Data.states{ts};
% ws     = state.wellSol;
% 
% isInj  = vertcat(ws.sign) > 0;
% ws_inj  = ws(isInj);
% ws_prod = ws(~isInj);
% % only selected
% ws_inj  = ws_inj(injIx);
% ws_prod = ws_prod(prodIx);
% 
% D = d.Data.diagnostics(ts).D;
% 
% if strcmp(mode, 'forward')
%     tri = sum(D.itracer(:,injIx), 2);
%     trp = D.ptracer(:,prodIx);
% elseif strcmp(mode, 'backward')
%     [ws_inj, ws_prod] = deal(ws_prod, ws_inj);
%     tri = sum(D.ptracer(:,prodIx), 2);
%     trp = D.itracer(:,injIx);
%     state.flux = -state.flux;
%     for k = 1:numel(ws_inj)
%         ws_inj(k).flux  = -ws_inj(k).flux;
%     end
%     for k = 1:numel(ws_prod)
%         ws_prod(k).flux = -ws_prod(k).flux;
%     end
% end        
% 
% % select subregion
% 
% tr = tri.*sum(trp, 2);
% cix = tr > 1e-5;
% if ~any(cix)
%     [t, qr, pvr, q, pv] = deal([]);
%     return;
% end
% [tri, trp]  =deal(tri(cix), trp(cix,:));
% tr = tr(cix);
% remap = zeros(d.G.cells.num, 1);
% remap(cix) = (1:nnz(cix))';
% remapplus = [0; remap];
% if isprop(d, 'Gs')
%     N = d.Gs.faces.neighbors;
% else
%     N = d.G.faces.neighbors;
% end
% N = remapplus(N+1);
% ie = prod(N,2) > 0;
% N = N(ie,:);
% v = state.flux(ie);
% 
% qi = sparse(vertcat(ws_inj.cells), 1, sum(vertcat(ws_inj.flux), 2), d.G.cells.num, 1);
% %qi = qi(cix).*tr;
% qi = qi(cix);
% ncp = arrayfun(@(x)numel(x.cells), ws_prod);
% qp = sparse(vertcat(ws_prod.cells), rldecode((1:numel(ws_prod))', ncp), sum(vertcat(ws_prod.flux), 2),  ...
%             d.G.cells.num, numel(ws_prod));
% qp = bsxfun(@times, qp(cix,:).*trp, tri);
% %qp = bsxfun(@times, qp(cix,:), tri);
% %qp = qp(cix,:);
% %% 
% s = zeros(nnz(cix),1);
% 
% upIx = zeros(size(v));
% upIx(v>=0) = N(v>=0, 1);
% upIx(v<0) = N(v<0, 2);
% 
% nf = size(N,1);
% nc = numel(s);
% C  = sparse( [(1:nf)'; (1:nf)'], N, ones(nf,1)*[1 -1], nf, nc);
% 
% pvfull  = d.Data.static(strcmp({d.Data.static.name}, 'PORV')).values;
% pvr  = pvfull(cix).*tr;
% %pvr  = pvfull(cix);
% 
% % select 200 time-steps from 0 tp 1*pvi
% nSteps = 100;
% %dt = 2*sum(pvr.*tr)/(nSteps*(sum(qi.*tr)+sum(sum(abs(qp.*trp),2))));
% dt = 2*sum(pvr.*tr)/(nSteps*(sum(qi.*tr)+sum(sum(abs(qp),2))));
% 
% s = qi./pvr;
% s = initVariablesADI(s);
% 
% eq = (C'*(v.*s(upIx).*tr(upIx))-s.*sum(qp,2))./pvr;
% %eq = (C'*(v.*s(upIx))-s.*sum(qp,2))./pvr;
% %eq = (C'*(v.*s(upIx))-s.*sum(qp,2))./pvr;
% s = double(s);
% A  = - eq.jac{1};
% M = (speye(numel(s)) -dt*A);
% q = zeros(numel(ws_prod), nSteps);
% for k = 1:nSteps
%     s = M\s;
%     q(:,k) = qp'*s;
% end
% q = q*year;
% q = bsxfun(@times, q, 1./sum(qp)');
% % replace 0/0 by 0
% q(~isfinite(q)) = 0;
% qr = sum(q,1);
% 
% %qr = qr/sum(sum(qp));
% pvr = sum(pvr)/(sum(qi)*year);
% pv  = bsxfun(@times, pvfull(cix).*tri, trp);
% pv  = sum(pv)/(sum(qi)*year);
% t = (1:nSteps)'*dt/year;
% 
% end
% 
% % 
% % 
% % dv = C'*v;
% % expr = @(sn, qi, qp)sn+dt*(C'*(v.*sn(upIx))-dv.*sn)./pvr + dt*(qi + sn.*sum(qp,2));
% % q = zeros(numel(ws_prod), nSteps);
% % for k = 1:nSteps
% %     s = expr(s, qi, qp);
% %     rem = s-1;
% %     while any(rem > sqrt(eps)) 
% %         s = expr(s, rem, 0);
% %         rem = s-1;
% %     end
% %     q(:,k) = qp'*s;
% % end
% % 
% % q = diff(q, 2)/dt;
% % qr = sum(q);
% % q = bsxfun(@mldivide, q, sum(qp)');
% % qr = qr/sum(sum(qp));
% % 
% % pvr = sum(pvr);
% % pv  = bsxfun(@mtimes, pvfull(cix).*sum(D.itracer(cix,injIx), 2), D.ptracer(cix,prodIx));
% % 
% % t = (1:nSteps-1)*dt;
% % end
