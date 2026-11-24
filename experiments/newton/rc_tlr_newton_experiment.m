function [results, firstPass] = rc_tlr_newton_experiment(dataDir, nValues, pFull, pLong)
%RC_TLR_NEWTON_EXPERIMENT RC-TLR on full and long-range Newton CP tensors.
%   p counts first-kind Chebyshev nodes (polynomial degree p-1).
%   Validation is relative maximum error on the complete middle slice.
%   Storage counts the selected core and all three evaluated grid factors.

    target = 1e-7;
    rrqrTol = 1e-9;
    strength = 2;

    if nargin < 1 || isempty(dataDir)
        dataDir = fullfile(fileparts(mfilename('fullpath')), ...
            'data','Newton_CP_-1_1');
    end
    if nargin < 2 || isempty(nValues)
        nValues = [256 512 1024 2048 4096];
    end
    if nargin < 3 || isempty(pFull)
        pFull = [129 257 513 1025 2049 4097 8193 16385];
    end
    if nargin < 4 || isempty(pLong)
        pLong = [129 257 513 1025 2049 4097];
    end
    validateattributes(nValues, {'numeric'}, ...
        {'vector','integer','>=',2,'finite','nonempty'});
    validateattributes(pFull, {'numeric'}, ...
        {'vector','integer','>=',2,'finite','nonempty'});
    validateattributes(pLong, {'numeric'}, ...
        {'vector','integer','>=',2,'finite','nonempty'});
    nValues = unique(nValues(:).');
    pFull = unique(pFull(:).');
    pLong = unique(pLong(:).');

    partNames = ["Full","Long"];
    pLists = {pFull,pLong};
    rows = cell(numel(nValues)*(numel(pFull)+numel(pLong)),13);
    summary = cell(2*numel(nValues),3);
    row = 0;
    group = 0;

    for n = nValues
        file = fullfile(dataDir,sprintf('Newt_canon_%d.mat',n));
        S = load(file,'U1r');
        assert(isfield(S,'U1r') && size(S.U1r,1)==n, ...
            'Expected an n-by-R array U1r in %s.',file);
        h = 2/n;
        x = -1 + ((1:n).'-0.5)*h;
        fullFactors = S.U1r/h;
        fullRank = size(fullFactors,2);
        longRank = ceil(fullRank/2)-3;
        assert(longRank>=1 && longRank<fullRank, ...
            'Invalid long-range split for n=%d.',n);

        for part = 1:2
            if part == 1
                P = fullFactors;
            else
                P = fullFactors(:,1:longRank);
            end
            firstP = NaN;
            for p = pLists{part}
                [C,F] = buildRC(P,x,p,rrqrTol,strength);
                err = sliceError(P,C,F);
                met = err <= target;
                if met && isnan(firstP), firstP = p; end
                r = size(F,2); % Symmetric input gives equal mode ranks.
                row = row+1;
                rows(row,:) = {partNames(part),n,p,size(P,2), ...
                    r,r,r,sprintf('%d x %d x %d',r,r,r),r, ...
                    numel(C),err,numel(C)+3*numel(F),met};
            end
            group = group+1;
            summary(group,:) = {partNames(part),n,firstP};
        end
    end
    results = cell2table(rows,'VariableNames', ...
        {'Part','n','p','InputCPRank','rX','rY','rZ','CoreDimensions', ...
         'TuckerRank','CoreStorage','RelativeError','TotalStorage','ToleranceMet'});
    firstPass = cell2table(summary,'VariableNames', ...
        {'Part','n','FirstPassingP'});
end

function [C, F] = buildRC(P, x, p, tol, strength)
    theta = ((2*(1:p).'-1)*pi)/(2*p);
    nodes = cos(theta);
    sqrtW = sqrt((pi/p)*sin(theta));
    samples = interp1(x,P,nodes,'spline','extrap');
    Z = sqrtW .* samples;

    % Represent the weighted unfolding in orthonormal row coordinates.
    % All singular directions are retained here. RRQR selects the rank.
    [~,S,V] = svd(Z,'econ');
    B = S*V.';
    q = size(B,1);

    % Symmetric packing preserves the Frobenius inner product.
    [i,j] = find(tril(true(q)));
    scale = ones(numel(i),1);
    scale(i~=j) = sqrt(2);
    pair = scale .* (B(i,:).*B(j,:));

    % Each column corresponds to a node in the selected mode.
    [I,J,T] = strongSelect(pair*Z.',tol,strength);
    r = numel(I);
    E = zeros(p,r);
    E(I,:) = eye(r);
    E(J,:) = T.' .* (sqrtW(I).' ./ sqrtW(J));

    C = symmetricCore(samples(I,:));
    baryWeights = (-1).^(0:p-1).' .* sin(theta);
    F = barycentricEval(nodes,baryWeights,E,x);
end

function [I, J, T] = strongSelect(A, tol, strength)
    [~,R,piv] = qr(A,0);
    if ~isvector(piv)
        [~,piv] = max(piv,[],1);
    end
    piv = piv(:);
    p = size(A,2);
    normA = norm(R,'fro');

    if normA == 0
        I = 1;
        J = (2:p).';
        T = zeros(1,p-1);
        return
    end

    R = R/normA;
    energy = sum(abs(R).^2,2);
    tails = [flipud(cumsum(flipud(energy(2:end)))); 0];
    r = find(tails <= tol^2,1,'first');
    swaps = 0;

    while true
        if r == p
            T = zeros(r,0);
            break
        end

        R11 = R(1:r,1:r);
        R22 = R(r+1:end,r+1:end);
        solved = R11 \ [R(1:r,r+1:end), eye(r)];
        T = solved(:,1:p-r);
        invR11 = solved(:,p-r+1:end);

        scores = abs(T).^2 + ...
            sum(abs(invR11).^2,2)*sum(abs(R22).^2,1);
        [maxScore, where] = max(scores(:));
        assert(isfinite(maxScore), 'Nonfinite strong RRQR score.');

        if maxScore > (strength*(1+1e-12))^2
            assert(swaps < 10*p, 'Strong RRQR exceeded its swap limit.');
            [i,j] = ind2sub(size(scores),where);
            j = r+j;
            R(:,[i j]) = R(:,[j i]);
            piv([i j]) = piv([j i]);
            [~,R] = qr(R,0);
            swaps = swaps+1;
        elseif norm(R22,'fro') <= tol
            break
        else
            assert(r < min(size(R)), 'RRQR could not meet its tolerance.');
            r = r+1;
        end
    end

    I = piv(1:r);
    J = piv(r+1:end);
end

function C = symmetricCore(A)
    % Sum the outer products of the CP columns.
    r = size(A,1);
    C = zeros(r,r,r);
    for k = 1:r
        C(:,:,k) = (A .* A(k,:))*A.';
    end
end

function F = barycentricEval(nodes, weights, values, x)
    % Apply the interpolation matrix in blocks to bound working memory.
    F = zeros(numel(x),size(values,2));
    block = 64;
    for first = 1:block:numel(x)
        ix = first:min(first+block-1,numel(x));
        D = x(ix)-nodes.';
        [distance, nearest] = min(abs(D),[],2);
        hit = distance == 0;

        local = zeros(numel(ix),size(values,2));
        local(hit,:) = values(nearest(hit),:);
        use = ~hit;
        if any(use)
            W = (distance(use)./D(use,:)).*weights.';
            local(use,:) = (W*values)./sum(W,2);
        end
        F(ix,:) = local;
    end
end

function err = sliceError(P,C,F)
    n = size(P,1);
    k = floor(n/2);
    r = size(F,2);
    sliceCore = reshape(reshape(C,r^2,r)*F(k,:).',r,r);
    maxRef = 0;
    maxErr = 0;
    for first = 1:128:n
        ix = first:min(first+127,n);
        reference = (P(ix,:).*P(k,:))*P.';
        approx = (F(ix,:)*sliceCore)*F.';
        maxRef = max(maxRef,max(abs(reference),[],'all'));
        maxErr = max(maxErr,max(abs(approx-reference),[],'all'));
    end
    assert(maxRef>0,'The reference slice has zero maximum magnitude.');
    err = maxErr/maxRef;
end
