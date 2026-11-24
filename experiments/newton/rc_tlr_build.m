function [model, stats] = rc_tlr_build(CP, grids, domain, p, qrTol, strength)
%RC_TLR_BUILD Construct the selected core, three maps, and three grid bases.
%   The core contains the unweighted selected Chebyshev values. The matrices
%   U,V,W (n_d-by-p) and E1,E2,E3 (p-by-r_d) are retained separately.
    [Gw, xc, yc, zc, G] = weighted_cheb_tensor(p, CP, grids, domain);
    I = cell(1,3);
    E = cell(1,3);
    basis = cell(1,3);
    stats = cell(1,3);
    referenceNodes = my_chev2(p);
    theta = (2*(1:p)-1)*pi/(2*p);
    baryWeights = (-1).^(0:p-1).*sin(theta);

    for d = 1:3
        other = setdiff(1:3, d, 'stable');
        A = Gw.U{other(1)};
        B = Gw.U{other(2)};
        Rcp = numel(Gw.lambda);
        pair = zeros(size(A,1)*size(B,1), Rcp);
        for k = 1:Rcp
            pair(:,k) = reshape(A(:,k)*B(:,k).', [], 1);
        end
        coeff = Gw.lambda(:) .* Gw.U{d}.';
        [I{d}, Ew, stats{d}] = direct_rrqr(pair, coeff, qrTol, strength);
        s = sqrt(Gw.weights{d}(:));
        E{d} = Ew.' .* (s(I{d}).' ./ s);
        E{d}(I{d},:) = eye(numel(I{d}));
        clear pair
    end

    C = selected_core(p, I, G);
    for d = 1:3
        grid = grids{d}(:);
        if isequal(domain(d,:), [-1 1])
            x = grid;
        else
            x = 2*(grid-domain(d,1))/(domain(d,2)-domain(d,1))-1;
        end
        basis{d} = zeros(numel(x),p);
        for first = 1:256:numel(x)
            rows = first:min(first+255,numel(x));
            basis{d}(rows,:) = lagrange_cheb(referenceNodes,x(rows),baryWeights);
        end
    end
    model = struct('C', C, 'E1', E{1}, 'E2', E{2}, 'E3', E{3}, ...
        'U', basis{1}, 'V', basis{2}, 'W', basis{3}, ...
        'selectedIndices', {I}, 'nodes', {{xc,yc,zc}}, ...
        'grids', {grids}, 'domain', domain, 'p', p, 'ranks', cellfun(@numel,I));
end

function [I, E, info] = direct_rrqr(pair, coeff, tol, strength)
% Each reflector is constructed from a pivot column of the full unfolding
% pair*coeff. Apply the reflector to pair to update all unfolding columns at
% once. The kth step zeros the subdiagonal in pivot k. There is no QR
% preprocessing of the CP factors or reduction of the selected core.
    K = pair;
    [rows,p] = deal(size(K,1),size(coeff,2));
    piv = 1:p;
    maxRank = min(rows,p);
    tails = product_column_norms(K,coeff);
    normA = norm(tails);
    swaps = 0;
    if normA == 0
        I = 1;
        E = zeros(1,p);
        E(1) = 1;
        residual = 0;
        maxRho = 0;
    else
        r = 0;
        while r < maxRank && (r == 0 || norm(tails) > tol*normA)
            [~, j] = max(tails);
            j = r+j;
            r = r+1;
            piv([r j]) = piv([j r]);
            tails = qr_step(r);
        end

        while true
            residual = norm(tails)/normA;
            maxRho = 0;
            T = zeros(r,0);
            if r < p
                top = (K(1:r,:)/normA)*coeff(:,piv);
                R11 = triu(top(:,1:r));
                F = R11 \ [top(:,r+1:end), eye(r)];
                T = F(:,1:p-r);
                inverse = F(:,p-r+1:end);
                rho2 = abs(T).^2 + sum(abs(inverse).^2,2)*((tails/normA).^2);
                [rhoMax2, pos] = max(rho2(:));
                maxRho = sqrt(rhoMax2);
                assert(isfinite(maxRho),'rc_tlr_build:NonfiniteRRQR', ...
                    'Nonfinite strong RRQR test at rank %d.',r);
                if maxRho > strength*(1+1e-12)
                    assert(swaps < 10*p,'rc_tlr_build:SwapLimit','Strong RRQR swap limit reached.');
                    [i,j] = ind2sub(size(rho2),pos);
                    j = r+j;
                    piv([i j]) = piv([j i]);
                    swaps = swaps+1;
                    K = pair;
                    for k = 1:r
                        tails = qr_step(k);
                    end
                    continue
                end
            end
            if residual <= tol, break; end
            assert(r < maxRank,'rc_tlr_build:Residual','Could not reach the RRQR tolerance.');
            [~, j] = max(tails);
            j = r+j;
            r = r+1;
            piv([r j]) = piv([j r]);
            tails = qr_step(r);
        end
        I = piv(1:r).';
        E = zeros(r,p);
        E(:,I) = eye(r);
        E(:,piv(r+1:end)) = T;
    end
    info = struct('relativeResidual',residual,'maxRho',maxRho, ...
        'numSwaps',swaps,'unfoldingSize',[rows p]);

    function remaining = qr_step(k)
        column = K(k:end,:)*coeff(:,piv(k));
        len = norm(column);
        assert(len > 0 && isfinite(len),'rc_tlr_build:ZeroPivot','Invalid QR pivot at step %d.',k);
        v = column/len;
        if v(1) < 0, v(1) = v(1)-1; else, v(1) = v(1)+1; end
        v = v/norm(v);
        change = 2*(v.'*K(k:end,:));
        for first = k:16384:rows
            ix = first:min(first+16383,rows);
            K(ix,:) = K(ix,:) - v(ix-k+1)*change;
        end
        remaining = product_column_norms(K(k+1:end,:),coeff(:,piv(k+1:end)));
    end
end

function norms = product_column_norms(A,B)
% Use the Gram of the current trailing rows for column lengths only. Directly
% accumulate norms in blocks when cancellation can affect a squared norm.
    if isempty(B) || isempty(A)
        norms = zeros(1,size(B,2));
        return
    end
    gram = A.'*A;
    squares = sum(B.*(gram*B),1);
    scales = sqrt(max(diag(gram),0)).'*abs(B);
    suspect = find(scales > 0 & squares <= 1e-10*scales.^2);
    for first = 1:16:numel(suspect)
        j = suspect(first:min(first+15,numel(suspect)));
        squares(j) = 0;
        for row = 1:16384:size(A,1)
            block = A(row:min(row+16383,size(A,1)),:)*B(:,j);
            squares(j) = squares(j)+sum(abs(block).^2,1);
        end
    end
    norms = sqrt(max(squares,0));
end
