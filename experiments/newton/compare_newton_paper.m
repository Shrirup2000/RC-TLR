function [results, firstPass] = compare_newton_paper( ...
    dataDir, nValues, pFull, pLong, repeats)
%COMPARE_NEWTON_PAPER RC-TLR and ChebTuck on the paper's Newton datasets.
%   P is the number of Chebyshev nodes. The polynomial degree is P-1.
%   Grid storage counts the core and all three evaluated factor matrices.
%   Build time covers formation of the core and evaluated grid factors.
%   Validation uses the full middle slice, processed in row blocks.

    if nargin < 2 || isempty(nValues)
        nValues = [256 512 1024 2048 4096];
    end
    if nargin < 3 || isempty(pFull)
        pFull = [129 257 513 1025 2049 4097 8193 16385];
    end
    if nargin < 4 || isempty(pLong)
        pLong = [129 257 513 1025 2049 4097];
    end
    if nargin < 5 || isempty(repeats)
        repeats = 1;
    end

    validateattributes(nValues, {'numeric'}, ...
        {'vector','integer','positive','finite'});
    validateattributes(pFull, {'numeric'}, ...
        {'vector','integer','>=',2,'finite'});
    validateattributes(pLong, {'numeric'}, ...
        {'vector','integer','>=',2,'finite'});
    validateattributes(repeats, {'numeric'}, ...
        {'scalar','integer','positive','finite'});

    nValues = unique(nValues(:).');
    pFull = unique(pFull(:).');
    pLong = unique(pLong(:).');

    target = 1e-7;
    rrqrTol = 1e-9;
    chebTol = 1e-7;
    strength = 2;

    partNames = ["Full", "Long"];
    pLists = {pFull, pLong};

    rows = cell(numel(nValues)*(numel(pFull)+numel(pLong)), 13);
    summary = cell(2*numel(nValues), 4);
    row = 0;
    group = 0;

    for n = nValues
        file = fullfile(dataDir, sprintf('Newt_canon_%d.mat', n));
        S = load(file, 'U1r');
        assert(isfield(S, 'U1r') && size(S.U1r,1) == n, ...
            'Expected an n-by-R array U1r in %s.', file);

        h = 2/n;
        x = (-1 + ((1:n).'-0.5)*h);
        fullFactors = S.U1r/h;
        fullRank = size(fullFactors,2);
        longRank = ceil(fullRank/2)-3;
        assert(longRank >= 1 && longRank < fullRank, ...
            'Invalid long-range split for n=%d.', n);

        for part = 1:2
            if part == 1
                P = fullFactors;
            else
                P = fullFactors(:,1:longRank);
            end

            firstRC = NaN;
            firstCT = NaN;

            for p = pLists{part}
                [Crc, Frc, timeRC] = timedBuild( ...
                    @() buildRC(P,x,p,rrqrTol,strength), repeats);
                [Cct, Fct, timeCT] = timedBuild( ...
                    @() buildCT(P,x,p,chebTol), repeats);

                [errRC, errCT] = sliceErrors(P,Crc,Frc,Cct,Fct);
                passRC = errRC <= target;
                passCT = errCT <= target;
                if passRC && isnan(firstRC)
                    firstRC = p;
                end
                if passCT && isnan(firstCT)
                    firstCT = p;
                end

                % Symmetric Newton data give equal ranks in the three modes.
                rankRC = size(Frc,2);
                rankCT = size(Fct,2);
                storageRC = numel(Crc) + 3*numel(Frc);
                storageCT = numel(Cct) + 3*numel(Fct);

                row = row+1;
                rows(row,:) = {partNames(part), n, p, size(P,2), ...
                    errRC, errCT, passRC, passCT, rankRC, rankCT, ...
                    storageRC, storageCT, [timeRC timeCT]};
            end

            group = group+1;
            summary(group,:) = {partNames(part), n, firstRC, firstCT};
        end
    end

    results = cell2table(rows(:,1:12), 'VariableNames', ...
        {'Part','n','p','InputCPRank','RC_Error','CT_Error', ...
         'RC_Pass','CT_Pass','RC_Rank','CT_Rank','RC_Storage','CT_Storage'});
    times = vertcat(rows{:,13});
    results.RC_BuildSeconds = times(:,1);
    results.CT_BuildSeconds = times(:,2);

    firstPass = cell2table(summary, 'VariableNames', ...
        {'Part','n','RC_FirstPassingP','CT_FirstPassingP'});
end

function [C, F, seconds] = timedBuild(builder, repeats)
    % Warm up this method and case before recording timings.
    builder();
    times = zeros(repeats,1);
    for k = 1:repeats
        start = tic;
        [C,F] = builder();
        times(k) = toc(start);
    end
    seconds = median(times);
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

function [C, F] = buildCT(P, x, p, tol)
    nodes = cos((0:p-1).'*pi/(p-1));
    samples = interp1(x,P,nodes,'spline','extrap');

    % DCT-I for Chebyshev coefficients via an even FFT extension.
    extended = [samples; samples(end-1:-1:2,:)];
    coeff = real(fft(extended,[],1))/(p-1);
    coeff = coeff(1:p,:);
    coeff([1 end],:) = coeff([1 end],:)/2;

    % Authors' RHOSVD rank rule for the three-mode CP coefficient tensor.
    [U,S,V] = svd(coeff,'econ');
    sigma = diag(S);
    tails = sqrt(flipud(cumsum(flipud(sigma.^2))));
    r = find(tails > tol/3,1,'last');
    assert(~isempty(r), 'ChebTuck found zero rank at this tolerance.');

    projected = sigma(1:r) .* V(:,1:r).';
    C = symmetricCore(projected);

    % Recover nodal values of the retained polynomial factors.
    kept = U(:,1:r);
    kept([1 end],:) = 2*kept([1 end],:);
    extended = [kept; kept(end-1:-1:2,:)];
    values = real(ifft((p-1)*extended,[],1));

    baryWeights = (-1).^(0:p-1).';
    baryWeights([1 end]) = baryWeights([1 end])/2;
    F = barycentricEval(nodes,baryWeights,values(1:p,:),x);
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

function [errRC, errCT] = sliceErrors(P, Crc, Frc, Cct, Fct)
    n = size(P,1);
    k = floor(n/2);
    rRC = size(Frc,2);
    rCT = size(Fct,2);

    sliceRC = reshape(reshape(Crc,rRC^2,rRC)*Frc(k,:).',rRC,rRC);
    sliceCT = reshape(reshape(Cct,rCT^2,rCT)*Fct(k,:).',rCT,rCT);

    maxRef = 0;
    maxRC = 0;
    maxCT = 0;
    for first = 1:128:n
        ix = first:min(first+127,n);
        reference = (P(ix,:).*P(k,:))*P.';
        approxRC = (Frc(ix,:)*sliceRC)*Frc.';
        approxCT = (Fct(ix,:)*sliceCT)*Fct.';

        maxRef = max(maxRef,max(abs(reference),[],'all'));
        maxRC = max(maxRC,max(abs(approxRC-reference),[],'all'));
        maxCT = max(maxCT,max(abs(approxCT-reference),[],'all'));
    end

    assert(maxRef > 0, 'The reference slice has zero maximum magnitude.');
    errRC = maxRC/maxRef;
    errCT = maxCT/maxRef;
end