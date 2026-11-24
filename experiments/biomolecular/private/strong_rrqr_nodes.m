function [I,selectedNodes,r,info] = strong_rrqr_nodes(Gw,mode,tol,f,previous,requestedRank)
%STRONG_RRQR_NODES Strong RRQR of a dense or CP weighted mode unfolding.
%   CP path uses an exact orthogonal reduction, not a Gram matrix.
%   No p-by-p-by-p sample array or p^2-by-p unfolding is formed on this path.
%   Supports rectangular reduced matrices and optional prescribed rank.
%   PREVIOUS may only be reused when its sourceMatrix equals the current one.
    if nargin < 3 || isempty(tol), tol = 1e-8; end
    if nargin < 4 || isempty(f), f = 2; end
    if nargin < 5, previous = []; end
    if nargin < 6, requestedRank = []; end
    validateattributes(mode,{'numeric'},{'scalar','integer','>=',1,'<=',3});
    validateattributes(tol,{'numeric'},{'scalar','real','finite','>',0,'<',1});
    validateattributes(f,{'numeric'},{'scalar','real','finite','>',1});
    if isstruct(Gw)
        p = size(Gw.U{mode},1);
        other = setdiff(1:3,mode,'stable');
        % Z2 (x) Z1 = (Q2 kron Q1)*(R2 (x) R1), columnwise product.
        [~,R1] = qr(Gw.U{other(1)},0);
        [~,R2] = qr(Gw.U{other(2)},0);
        Rcp = numel(Gw.lambda);
        pair = reshape(reshape(R1,size(R1,1),1,Rcp) .* ...
            reshape(R2,1,size(R2,1),Rcp),[],Rcp);
        [~,Rpair] = qr(pair,0);
        A = (Rpair .* Gw.lambda(:).') * Gw.U{mode}.';
        if isfield(Gw,'nodes'), nodes = Gw.nodes{mode};
        else, nodes = my_chev2(p); end
        if isfield(Gw,'weights'), weights = Gw.weights{mode};
        else, weights = []; end
        route = 'CP orthogonal reduction';
    else
        dims = [size(Gw,1),size(Gw,2),size(Gw,3)]; p = dims(mode);
        order = [mode,setdiff(1:3,mode,'stable')];
        A = reshape(permute(Gw,order),p,[]).';
        nodes = my_chev2(p); weights = []; route = 'dense unfolding';
    end
    assert(~isempty(A) && isreal(A) && all(isfinite(A(:))), ...
        'The unfolding must be nonempty, real and finite.');
    maxRank = min(size(A));
    if ~isempty(requestedRank)
        validateattributes(requestedRank,{'numeric'}, ...
            {'scalar','integer','positive','<=',maxRank});
    end
    if isempty(previous)
        [~,R,P] = qr(A,0);
        if isvector(P), piv = P(:).'; else, [~,piv] = max(P,[],1); end
    else
        assert(isfield(previous,'sourceMatrix') && previous.mode==mode && ...
            isequal(previous.sourceMatrix,A), ...
            'Previous factorization is not for the same tensor and mode.');
        R = previous.R; piv = previous.piv;
    end
    normA = norm(R,'fro'); numSwaps = 0; maxSwaps = 10*p;
    strongLimit = f*(1+1e-12);
    if normA == 0
        r = 1; E = zeros(1,p); E(piv(1)) = 1;
        relResidual = 0; maxRho = 0;
    else
        % Scaling improves the triangular solves; it does not change E/rho.
        R = R/normA;
        if isempty(requestedRank)
            energy = sum(abs(R).^2,2);
            tails = [flipud(cumsum(flipud(energy(2:end))));0];
            r = find(tails<=tol^2,1,'first');
        else
            r = requestedRank;
        end
        while true
            if r == p
                maxRho = 0; relResidual = 0; break
            end
            R11 = R(1:r,1:r);
            if rcond(R11)==0
                error('strong_rrqr_nodes:SingularBlock', ...
                    'Selected block is singular. Reduce prescribed rank or relax QR tolerance.');
            end
            B = R11 \ [R(1:r,r+1:p),eye(r)];
            T = B(:,1:p-r); invR11 = B(:,p-r+1:end);
            R22 = R(r+1:end,r+1:p);
            rho2 = abs(T).^2 + sum(abs(invR11).^2,2)*sum(abs(R22).^2,1);
            [rhoMax2,index] = max(rho2(:)); maxRho = sqrt(rhoMax2);
            if ~isfinite(maxRho)
                error('strong_rrqr_nodes:Nonfinite','Numerically singular strong-RRQR test.');
            end
            if maxRho > strongLimit
                if numSwaps >= maxSwaps
                    error('strong_rrqr_nodes:SwapLimit','Strong RRQR exceeded %d swaps.',maxSwaps);
                end
                [i,j] = ind2sub(size(rho2),index); j = j+r;
                R(:,[i j]) = R(:,[j i]); piv([i j]) = piv([j i]);
                [~,R] = qr(R,0);
                numSwaps = numSwaps+1;
                continue
            end
            relResidual = norm(R22,'fro');
            if ~isempty(requestedRank) || relResidual <= tol, break; end
            if r==maxRank
                error('strong_rrqr_nodes:Residual','Cannot meet tolerance at available rank.');
            end
            r = r+1;
        end
        E = zeros(r,p); E(:,piv(1:r)) = eye(r);
        if r<p, E(:,piv(r+1:p)) = R(1:r,1:r) \ R(1:r,r+1:p); end
        R = R*normA;
    end
    I = piv(1:r).'; selectedNodes = reshape(nodes(I),[],1);
    info = struct('R',R,'piv',piv,'rank',r,'normA',normA,'mode',mode, ...
        'p',p,'relativeResidual',relResidual,'maxRho',maxRho, ...
        'numSwaps',numSwaps,'accuracySatisfied',relResidual<=tol, ...
        'strongSatisfied',maxRho<=strongLimit,'E',E,'weights',weights, ...
        'sourceMatrix',A,'route',route,'reducedSize',size(A));
end
