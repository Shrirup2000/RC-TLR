function [results, comparison] = rc_tlr_threebox_experiment( ...
    p_values, tol_values, n, f, pointSpec, coreMethod, stopOnAll, rrqrTol)
%RC_TLR_THREEBOX_EXPERIMENT Weighted and unweighted RC-TLR on three boxes.
%   Default p=2:8, target 1e-7, uniform n=20 per axis, strength f=2,
%   strong RRQR tolerance 4e-8. The full fixed sweep is run by default.
%   Both variants use first-kind nodes and the unweighted selected core.
%   TotalTime is construction only (sampling, selected core, maps, 1D bases).
%   Target reference evaluation and reconstruction are excluded from timing.
%   TotalStorage = rX*rY*rZ + p^2*(rX+rY+rZ) + 6*n*p for Cartesian grids.
%   Pass [] for rrqrTol to use Tol/sqrt(3), or true for stopOnAll to stop
%   when both variants meet the target. All helpers are local functions.

    if nargin < 1 || isempty(p_values), p_values = 2:8; end
    validateattributes(p_values, {'numeric'}, ...
        {'vector','real','finite','integer','>=',2}, mfilename, 'p_values');
    p_values = double(p_values(:));
    if nargin < 2 || isempty(tol_values), tol_values = 1e-7; end
    if nargin < 3 || isempty(n), n = 20; end
    if nargin < 4 || isempty(f), f = 2; end
    if nargin < 5 || isempty(pointSpec), pointSpec = 'uniform'; end
    if nargin < 6 || isempty(coreMethod), coreMethod = 'skeleton'; end
    if nargin < 7 || isempty(stopOnAll), stopOnAll = false; end
    if nargin < 8, rrqrTol = 4e-8; end

    validateattributes(tol_values, {'numeric'}, ...
        {'vector','real','finite','>',0,'<',1}, mfilename, 'tol_values');
    validateattributes(n, {'numeric'}, {'scalar','finite','integer','>=',2}, ...
        mfilename, 'n');
    validateattributes(f, {'numeric'}, {'scalar','real','finite','>',1}, ...
        mfilename, 'f');
    n = double(n);
    f = double(f);
    validateattributes(stopOnAll, {'logical','numeric'}, {'scalar','binary'}, ...
        mfilename, 'stopOnAll');
    assert(strcmpi(coreMethod,'skeleton'), 'coreMethod must be ''skeleton''.');
    if ~isempty(rrqrTol)
        validateattributes(rrqrTol, {'numeric'}, ...
            {'scalar','real','finite','>',0,'<',1}, mfilename, 'rrqrTol');
    end

    tol_values = double(tol_values(:));
    if isscalar(p_values)
        p_values = repmat(p_values, numel(tol_values), 1);
    elseif isscalar(tol_values)
        tol_values = repmat(tol_values, numel(p_values), 1);
    end
    assert(numel(p_values) == numel(tol_values), ...
        'Supply one p per tolerance, or make either input scalar.');

    boxes = {[-3,-1;0,2], [1,3;0,2], [-1,1;4,6]};
    if isstruct(pointSpec)
        assert(isscalar(pointSpec) && all(isfield(pointSpec, {'X','Y','Z'})), ...
            'pointSpec must have X, Y and Z fields.');
        points = {pointSpec.X, pointSpec.Y, pointSpec.Z};
        isCartesian = false;
    else
        distribution = lower(char(pointSpec));
        points = cell(1,3);
        [points{1:3}] = uniform_points(n,distribution);
        isCartesian = ismember(distribution, {'uniform','grid','random_tensor'});
    end

    basisPoints = cell(3,2);
    for mode = 1:3
        validateattributes(points{mode}, {'double','single'}, ...
            {'2d','ncols',2,'nonempty','real','finite'}, mfilename, 'particles');
        points{mode} = double(points{mode});
        box = boxes{mode};
        margin = 16*eps(max(1,max(abs(box(:)))));
        assert(all(all(points{mode} >= box(:,1).'-margin & ...
            points{mode} <= box(:,2).'+margin)), ...
            'Particle coordinates must lie inside their box.');
        if isCartesian
            x1 = reshape(points{mode}(:,1), n, n);
            x2 = reshape(points{mode}(:,2), n, n);
            basisPoints(mode,:) = {x1(:,1), x2(1,:).'};
        else
            basisPoints(mode,:) = {points{mode}(:,1), points{mode}(:,2)};
        end
    end

    Dxy = squared_dist(points{1},points{2});
    Dyz = squared_dist(points{2},points{3});
    Dxz = squared_dist(points{1},points{3});
    Nz = size(points{3},1);
    blockSize = 16;
    normK2 = 0;
    for k0 = 1:blockSize:Nz
        ids = k0:min(k0+blockSize-1,Nz);
        Ktrue = exact_kernel_block(Dxy,Dyz,Dxz,ids);
        normK2 = normK2 + sum(Ktrue(:).^2);
    end
    clear Ktrue

    names = ["Weighted","Unweighted"];
    rowTemplate = struct('Tol',0,'p',0,'Variant',"", ...
        'rX',0,'rY',0,'rZ',0,'RankSum',0,'CoreStorage',0, ...
        'TotalStorage',0,'TotalTime',0,'RelativeError',0, ...
        'ToleranceMet',false,'CompressionTolerance',0);
    rows = repmat(rowTemplate,2*numel(p_values),1);
    nodeOverlap = zeros(numel(p_values),3);

    fprintf('\nRC-TLR (weighted/unweighted)\n');
    fprintf('Particles per box: %d, %d, %d. TotalTime is construction only.\n', ...
        size(points{1},1),size(points{2},1),Nz);

    stoppedEarly = false;
    for ip = 1:numel(p_values)
        p = p_values(ip);
        tol = tol_values(ip);
        modeTol = rrqrTol;
        if isempty(modeTol), modeTol = tol/sqrt(3); end
        % Sample the shared RC core once and charge that cost to both variants.
        tSample = tic;
        G = chebyshev_core(p);
        sampleTime = toc(tSample);

        for iv = 1:2
            tBuild = tic;
            [C,E,selected] = rc_skeleton(G,p,modeTol,f,iv == 1);
            B = evaluation_bases(basisPoints,boxes,p);
            compressionTolerance = modeTol;
            totalTime = toc(tBuild)+sampleTime;

            ranks = cellfun(@(factor) size(factor,2), E);
            totalStorage = numel(C) + sum(cellfun(@numel,E)) + ...
                sum(cellfun(@numel,B(:)));

            % Match the selected node sets before leaving the two RC runs.
            if iv == 1
                weightedNodes = selected;
            elseif iv == 2
                for mode = 1:3
                    nodeOverlap(ip,mode) = numel(intersect( ...
                        weightedNodes{mode},selected{mode}));
                end
                clear G
            end

            relativeError = particle_error(C,E,B,isCartesian, ...
                Dxy,Dyz,Dxz,normK2,blockSize);
            row = 2*(ip-1) + iv;
            rows(row) = struct('Tol',tol,'p',p,'Variant',names(iv), ...
                'rX',ranks(1),'rY',ranks(2),'rZ',ranks(3), ...
                'RankSum',sum(ranks),'CoreStorage',numel(C), ...
                'TotalStorage',totalStorage,'TotalTime',totalTime, ...
                'RelativeError',relativeError,'ToleranceMet',relativeError <= tol, ...
                'CompressionTolerance',compressionTolerance);
            fprintf(['p=%2d tol=%.1e %-10s core=%3d x %3d x %3d  ' ...
                'error=%.6e  storage=%d  construction=%.6f s\n'], ...
                p,tol,char(names(iv)),ranks(1),ranks(2),ranks(3), ...
                relativeError,totalStorage,totalTime);
            clear C E B
        end

        if stopOnAll && all([rows(row-1:row).ToleranceMet])
            fprintf('Both RC-TLR variants meet tolerance %.1e at p=%d.\n',tol,p);
            stoppedEarly = true;
            break;
        end
    end

    if stopOnAll && ~stoppedEarly
        fprintf('Sweep finished before both RC-TLR variants met the target at the same p.\n');
    end

    results = struct2table(rows(1:row));
    results.Variant = categorical(results.Variant);
    Rw = results(1:2:end,:);
    Ru = results(2:2:end,:);
    nodeOverlap = nodeOverlap(1:height(Rw),:);
    comparison = table(Rw.Tol,Rw.p, ...
        Rw.rX,Ru.rX,Rw.rY,Ru.rY,Rw.rZ,Ru.rZ, ...
        Rw.RelativeError,Ru.RelativeError,Ru.RelativeError./Rw.RelativeError, ...
        Rw.TotalStorage,Ru.TotalStorage,Rw.TotalTime,Ru.TotalTime, ...
        Rw.ToleranceMet,Ru.ToleranceMet, ...
        nodeOverlap(:,1),nodeOverlap(:,2),nodeOverlap(:,3), ...
        'VariableNames',{'Tol','p','rX_w','rX_u','rY_w','rY_u', ...
        'rZ_w','rZ_u','Err_w','Err_u','ErrRatio_u_over_w', ...
        'Storage_w','Storage_u','Time_w','Time_u','Met_w','Met_u', ...
        'sameNodesX','sameNodesY','sameNodesZ'});
    fprintf('\nSide-by-side results (weighted and unweighted RC-TLR)\n');
    disp(comparison);
    plot_variants(Rw,Ru);
end

function [C,E,selected] = rc_skeleton(G,p,modeTol,f,useWeights)
    if useWeights
        [WX,WY,WZ] = cheb_weight_matrices(p);
        weights = {WX,WY,WZ};
        Gin = weighted_core(G,WX,WY,WZ);
    else
        weights = repmat({ones(p*p,1)},1,3);
        Gin = G;
    end
    E = cell(1,3);
    selected = cell(1,3);
    for mode = 1:3
        [selected{mode},~,rank,info] = strong_rrqr_nodes(Gin,mode,modeTol,f);
        E{mode} = rrqr_interpolation(info.R,info.piv,rank,weights{mode});
    end
    C = G(selected{1},selected{2},selected{3});
end

function B = evaluation_bases(basisPoints,boxes,p)
    B = cell(3,2);
    for mode = 1:3
        for axis = 1:2
            interval = boxes{mode}(axis,:);
            B{mode,axis} = lagrange_cheb(interval(1),interval(2),p,basisPoints{mode,axis});
        end
    end
end

function F = evaluate_particle_factor(B1,B2,E,isCartesian)
    % Cartesian: F = kron(B2,B1)*E. Scattered: use each particle's basis row.
    p = size(B1,2);
    n1 = size(B1,1);
    rank = size(E,2);
    if isCartesian
        n2 = size(B2,1);
        T = reshape(B1*reshape(E,p,[]),n1,p,rank);
        F = reshape(mode_product(T,B2,2,3),n1*n2,rank);
    else
        F = zeros(n1,rank,'like',E);
        for b = 1:p
            ids = (b-1)*p + (1:p);
            F = F + (B1*E(ids,:)).*B2(:,b);
        end
    end
end

function relErr = particle_error(C,E,B,isCartesian,Dxy,Dyz,Dxz,normK2,bs)
    FX = evaluate_particle_factor(B{1,1},B{1,2},E{1},isCartesian);
    FY = evaluate_particle_factor(B{2,1},B{2,2},E{2},isCartesian);
    FZ = evaluate_particle_factor(B{3,1},B{3,2},E{3},isCartesian);
    HX = mode_product(C,FX,1,3);
    err2 = 0;
    for k0 = 1:bs:size(FZ,1)
        ids = k0:min(k0+bs-1,size(FZ,1));
        T = mode_product(HX,FZ(ids,:),3,3);
        Khat = mode_product(T,FY,2,3);
        Ktrue = exact_kernel_block(Dxy,Dyz,Dxz,ids);
        delta = Ktrue-Khat;
        err2 = err2 + sum(delta(:).^2);
    end
    relErr = sqrt(err2/normK2);
end

function D = squared_dist(A,B)
    D = (A(:,1)-B(:,1).').^2 + (A(:,2)-B(:,2).').^2;
end

function K = exact_kernel_block(Dxy,Dyz,Dxz,ids)
    [Nx,Ny] = size(Dxy);
    nb = numel(ids);
    K = 1./sqrt(Dxy + reshape(Dyz(:,ids),1,Ny,nb) + ...
        reshape(Dxz(:,ids),Nx,1,nb));
end

function T = mode_product(A,M,mode,d)
    dims = size(A,1:d);
    order = [mode,1:mode-1,mode+1:d];
    B = M*reshape(permute(A,order),dims(mode),[]);
    newDims = dims(order);
    newDims(1) = size(M,1);
    T = ipermute(reshape(B,newDims),order);
end

function plot_variants(Rw,Ru)
    if numel(unique(Rw.p)) > 1 || isscalar(Rw.p)
        x = Rw.p;
        xLabel = 'Chebyshev nodes per coordinate (p)';
        xScale = 'linear';
    else
        x = Rw.Tol;
        xLabel = 'Target relative error';
        xScale = 'log';
    end
    fields = {'RelativeError','TotalStorage','TotalTime'};
    labels = {'Relative Frobenius error','Scalars stored','Construction time (s)'};
    yScales = {'log','log','linear'};
    for k = 1:numel(fields)
        figure('Name',labels{k},'Color','w');
        ax = axes;
        field = fields{k};
        plot(ax,x,Rw.(field),'-o',x,Ru.(field),'-s','LineWidth',1.5);
        legendLabels = {'RC-TLR (weighted)','RC-TLR (unweighted)'};
        if k == 1
            hold(ax,'on');
            plot(ax,x,Rw.Tol,'k--');
            legendLabels{end+1} = 'Target';
        end
        set(ax,'XScale',xScale,'YScale',yScales{k});
        grid(ax,'on');
        xlabel(ax,xLabel);
        ylabel(ax,labels{k});
        legend(ax,legendLabels,'Location','best');
    end
end

function [X,Y,Z] = uniform_points(n,distribution)
%UNIFORM_POINTS Evaluate three boxes with selectable particle distributions.
%   [X,Y,Z] = uniform_points(n)            % Cartesian uniform grid
%   [X,Y,Z] = uniform_points(n,'uniform')  % same as above
%   [X,Y,Z] = uniform_points(n,'random')   % independent scattered points
%   [X,Y,Z] = uniform_points(n,'random_tensor')
%                                       % random tensor-product axes
%
% N=n^2 points are returned per box, in N-by-2 arrays.
% For repeatable randomness call rng(seed) before this function.
% Function name and first argument remain unchanged.

    if nargin<2 || isempty(distribution), distribution='uniform'; end
    validateattributes(n,{'numeric'}, ...
        {'scalar','integer','>=',2},mfilename,'n');
    bounds = cat(3,[-3,-1;0,2],[1,3;0,2],[-1,1;4,6]);
    points = cell(1,3);
    distribution = lower(char(distribution));
    for s=1:3
        box=bounds(:,:,s);
        switch distribution
            case {'uniform','grid'}
                u=linspace(box(1,1),box(1,2),n);
                v=linspace(box(2,1),box(2,2),n);
                [U,V]=ndgrid(u,v);
                pts=[U(:),V(:)];
            case {'random','scattered'}
                pts=[box(1,1)+(box(1,2)-box(1,1))*rand(n*n,1), ...
                     box(2,1)+(box(2,2)-box(2,1))*rand(n*n,1)];
            case 'random_tensor'
                u=sort(box(1,1)+(box(1,2)-box(1,1))*rand(n,1));
                v=sort(box(2,1)+(box(2,2)-box(2,1))*rand(n,1));
                [U,V]=ndgrid(u,v);
                pts=[U(:),V(:)];
            otherwise
                error('Unknown distribution: %s',distribution);
        end
        points{s}=pts;
    end
    [X,Y,Z]=deal(points{:});
end

function G = chebyshev_core(p)

    % Step 1: Chebyshev nodes for box X
    x1 = my_chev2(-3, -1, p);
    x2 = my_chev2( 0,  2, p);

    % Step 2: Chebyshev nodes for box Y
    y1 = my_chev2(1, 3, p);
    y2 = my_chev2(0, 2, p);

    % Step 3: Chebyshev nodes for box Z
    z1 = my_chev2(-1, 1, p);
    z2 = my_chev2( 4, 6, p);

    % Step 4: Construct tensor-product Chebyshev grids
    [X1, X2] = ndgrid(x1, x2);
    [Y1, Y2] = ndgrid(y1, y2);
    [Z1, Z2] = ndgrid(z1, z2);

    % Step 5: Convert grids into coordinate matrices
    X = [X1(:), X2(:)];
    Y = [Y1(:), Y2(:)];
    Z = [Z1(:), Z2(:)];

    % Step 6: Number of Chebyshev grid points per box
    q = p^2;

    % Step 7: Compute squared pairwise distances
    Dxy = (X(:,1)-Y(:,1)').^2 ...
        + (X(:,2)-Y(:,2)').^2;

    Dyz = (Y(:,1)-Z(:,1)').^2 ...
        + (Y(:,2)-Z(:,2)').^2;

    Dxz = (X(:,1)-Z(:,1)').^2 ...
        + (X(:,2)-Z(:,2)').^2;

    % Step 8: Preallocate core tensor
    G = zeros(q,q,q);

    % Step 9: Construct the core one slice at a time
    for k = 1:q

        G(:,:,k) = 1 ./ sqrt( ...
            Dxy + Dyz(:,k).' + Dxz(:,k));

    end

end

function [I, selectedNodes, r, info] = strong_rrqr_nodes( ...
    Gw, mode, tol, f, previous)

    if nargin < 4 || isempty(f)
        f = 2;
    end

    if nargin < 5
        previous = [];
    end

    m = size(Gw,1);
    p = round(sqrt(m));

    %% Step 1: QR factorization or reuse previous QR

    if isempty(previous)

        % Transposed mode unfolding
        order = [mode, setdiff(1:3,mode)];

        A = reshape(permute(Gw,order),m,[]).';

        % Initial column-pivoted QR
        [~,R,piv] = qr(A,'econ','vector');

        piv = piv(:).';

        % Original matrix norm
        normA = norm(R,'fro');

        % Minimum admissible rank
        rMin = 1;

    else

        % Reuse previous strong RRQR factorization
        R = previous.R;
        piv = previous.piv;
        normA = previous.normA;

        % Never decrease rank as tolerance becomes stricter
        rMin = previous.rank;

    end

    %% Step 2: Initial rank from the QR residual

    rowEnergy = sum(abs(R).^2,2);

    tailEnergy = [ ...
        flipud(cumsum(flipud(rowEnergy(2:end)))); 0];

    target2 = (tol*normA)^2;

    r0 = find(tailEnergy <= target2,1,'first');

    r = max(rMin,r0);

    %% Step 3: Strong RRQR refinement

    nSwaps = 0;
    maxSwaps = 10*m;

    while true

        if r == m
            relResidual = 0;
            maxRho = 0;
            break;
        end

        % Triangular blocks
        R11 = R(1:r,1:r);
        R12 = R(1:r,r+1:m);
        R22 = R(r+1:m,r+1:m);

        % Interpolation coefficients
        T = R11 \ R12;

        % Strong RRQR criterion
        S = R11 \ eye(r);

        rowNorm2 = sum(abs(S).^2,2);
        colNorm2 = sum(abs(R22).^2,1);

        rho2 = abs(T).^2 + rowNorm2*colNorm2;

        [rhoMax2,idx] = max(rho2(:));

        maxRho = sqrt(rhoMax2);

        %% Strong RRQR column swap

        if maxRho > f*(1+1e-12)

            if nSwaps >= maxSwaps
                error(['Strong RRQR refinement did not ' ...
                       'converge within the swap limit.']);
            end

            [i,j] = ind2sub(size(rho2),idx);

            j = r+j;

            % Exchange selected and unselected columns
            permutation = 1:m;
            permutation([i,j]) = permutation([j,i]);

            piv = piv(permutation);

            % Re-triangularize the small R factor
            R = qr(R(:,permutation));

            nSwaps = nSwaps+1;

            continue;
        end

        %% Step 4: Check residual tolerance

        relResidual = norm(R22,'fro')/normA;

        if relResidual <= tol
            break;
        end

        % Increase rank if tolerance is not satisfied
        r = r+1;

    end

    %% Step 5: Selected node indices

    I = piv(1:r).';

    %% Step 6: Selected Chebyshev coordinates

    switch mode

        case 1
            intervals = [-3,-1; 0,2];

        case 2
            intervals = [1,3; 0,2];

        case 3
            intervals = [-1,1; 4,6];

        otherwise
            error('Mode must be 1, 2 or 3.');

    end

    c1 = my_chev2(intervals(1,1),intervals(1,2),p);
    c2 = my_chev2(intervals(2,1),intervals(2,2),p);

    c1 = c1(:);
    c2 = c2(:);

    % ndgrid-compatible ordering
    i1 = mod(I-1,p)+1;
    i2 = floor((I-1)/p)+1;

    selectedNodes = [c1(i1),c2(i2)];

    %% Step 7: Save factorization for next tolerance

    info.R = R;
    info.piv = piv;
    info.rank = r;
    info.normA = normA;

    info.relativeResidual = relResidual;
    info.maxRho = maxRho;
    info.numSwaps = nSwaps;

    info.accuracySatisfied = relResidual <= tol;
    info.strongSatisfied = maxRho <= f*(1+1e-12);

end

function [P,L] = rrqr_interpolation(R,piv,r,W)

    % Number of Chebyshev nodes
    m = size(R,2);

    % Selected and unselected indices
    I = piv(1:r);
    J = piv(r+1:m);

    % Compute interpolation coefficients
    T = R(1:r,1:r) \ R(1:r,r+1:m);

    % Square roots of quadrature weights
    s = sqrt(W(:));

    % Allocate interpolation matrix
    P = zeros(m,r,'like',R);

    % Interpolation property at selected nodes
    P(I,:) = eye(r,'like',R);

    % Interpolation at unselected nodes
    P(J,:) = T.' .* (s(I).' ./ s(J));

    % Construct weighted interpolation matrix if requested
    if nargout > 1

        L = zeros(m,r,'like',R);

        L(I,:) = eye(r,'like',R);
        L(J,:) = T.';

    end

end

function [WX, WY, WZ, DX, DY, DZ] = cheb_weight_matrices(p)

    % Step 1: Generate Chebyshev indices
    k = (1:p)';

    % Step 2: Compute 1D Chebyshev integration weights
    w = (pi/p) * sin((2*k-1)*pi/(2*p));

    % Step 3: Tensor-product weights for 2D boxes
    W = kron(w, w);

    % Step 4: Weight vectors for X, Y and Z
    WX = W;
    WY = W;
    WZ = W;

    % Step 5: Construct sparse diagonal matrices only if requested
    if nargout > 3

        D = spdiags(W, 0, p^2, p^2);

        DX = D;
        DY = D;
        DZ = D;

    end

end

function L = lagrange_cheb(a, b, n, x)

    % Step 1: Generate Chebyshev nodes using my_chev2
    xi = my_chev2(a, b, n);

    % Step 2: Compute barycentric weights
    k = 1:n;
    theta = (2*k-1)*pi/(2*n);
    w = (-1).^(k-1).*sin(theta);

    % Step 3: Arrange evaluation points as a column
    x = x(:);
    m = numel(x);

    % Step 4: Preallocate Lagrange basis matrix
    L = zeros(m, n);

    % Step 5: Identify points coinciding with nodes
    [onNode, idx] = ismember(x, xi);

    % Step 6: Set exact interpolation values
    rows = find(onNode);

    if ~isempty(rows)
        L(sub2ind([m,n], rows, idx(onNode))) = 1;
    end

    % Step 7: Evaluate basis at remaining points
    rows = find(~onNode);

    if ~isempty(rows)

        D = x(rows) - xi;

        T = w ./ D;

        L(rows,:) = T ./ sum(T, 2);

    end

end

function Gw = weighted_core(G, WX, WY, WZ)

    % Step 1: Number of Chebyshev grid points
    m = size(G, 1);

    % Step 2: Compute square roots of weight vectors
    sx = sqrt(WX(:));
    sy = sqrt(WY(:));
    sz = sqrt(WZ(:));

    % Step 3: Combine weights of first two modes
    A = sx * sy.';

    % Step 4: Preallocate weighted core
    Gw = zeros(size(G), 'like', G);

    % Step 5: Apply weights slice by slice
    for k = 1:m

        Gw(:,:,k) = sz(k) * (A .* G(:,:,k));

    end

end

function X = my_chev2(a, b, n)

    % Step 1: Generate integer indices
    k = 1:n;

    % Step 2: Chebyshev nodes of the first kind on [-1,1]
    x_cheb = cos((2*k-1)*pi/(2*n));

    % Step 3: Map the nodes from [-1,1] to [a,b]
    X = (a+b)/2 + (b-a)/2*x_cheb;

end
