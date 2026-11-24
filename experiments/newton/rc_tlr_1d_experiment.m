function [results, model, diagnostics] = rc_tlr_1d_experiment(options)
%RC_TLR_1D_EXPERIMENT RC-TLR on the long-range Newton CP input.
%   The defaults are n=256, pValues=[129 257 513 1025], relative-error
%   target tol=1e-7, and weighted RRQR tolerance qrTol=5e-8. The default
%   validation metric is relative maximum error on the whole middle slice.
%   Set validation='full' for the full-volume relative Frobenius error.
%
%   TotalStorage counts the core and three evaluated grid factors.
%
%   Defaults: strength=2, blockSize=8, stopWhenMet=true.
%   For custom input, provide CP, grids={x,y,z}, and the physical domain.
    if nargin < 1 || isempty(options), options = struct; end
    assert(isstruct(options) && isscalar(options),'Options must be a scalar struct.');
    opts = struct('n',256,'pValues',[129 257 513 1025], ...
        'tol',1e-7,'qrTol',5e-8,'strength',2,'blockSize',8, ...
        'validation','slice','stopWhenMet',true, ...
        'dataDir',[],'CP',[],'grids',[],'domain',repmat([-1 1],3,1));
    names = fieldnames(options);
    for k = 1:numel(names)
        field = names{k};
        assert(isfield(opts,field),'Unknown option: %s',field);
        opts.(field) = options.(field);
    end
    validateattributes(opts.pValues,{'numeric'}, ...
        {'vector','integer','positive','finite','nonempty'});
    validateattributes(opts.tol,{'numeric'},{'scalar','real','finite','>',0,'<',1});
    validateattributes(opts.qrTol,{'numeric'},{'scalar','real','finite','>',0,'<',1});
    validateattributes(opts.strength,{'numeric'},{'scalar','real','finite','>',1});
    validateattributes(opts.blockSize,{'numeric'},{'scalar','integer','positive','finite'});
    validateattributes(opts.stopWhenMet,{'logical','numeric'},{'scalar','binary'});
    validateattributes(opts.domain,{'numeric'},{'size',[3 2],'real','finite'});
    assert(all(opts.domain(:,2)>opts.domain(:,1)),'Invalid domain.');
    opts.validation = validatestring(opts.validation,{'slice','full'});

    if isempty(opts.CP)
        assert(isequal(opts.domain,repmat([-1 1],3,1)), ...
            'The Newton data are defined on [-1,1]^3.');
        [CP,x,y,z] = kernel_tensor(opts.n,'cp','long',opts.dataDir);
        grids = {x,y,z};
    else
        CP = opts.CP;
        grids = opts.grids;
        assert(iscell(grids) && numel(grids)==3,'Custom CP input requires grids={x,y,z}.');
        grids = reshape(grids,1,3);
    end
    n = cellfun(@numel,grids);
    pValues = opts.pValues(:);
    count = numel(pValues);
    ranks = zeros(count,3);
    coreStorage = cell(count,1);
    relativeError = zeros(count,1);
    totalStorage = zeros(count,1);
    diagnostics = cell(count,1);

    fprintf('\nRC-TLR 1D experiment\n');
    fprintf('CP rank = %d, grid = %d x %d x %d, target = %.1e, RRQR tolerance = %.1e\n', ...
        numel(CP.lambda),n(1),n(2),n(3),opts.tol,opts.qrTol);
    if strcmp(opts.validation,'slice')
        fprintf('RelativeError uses the full middle slice and its maximum norm.\n');
    else
        fprintf('RelativeError uses the full-grid Frobenius norm.\n');
    end
    fprintf('%5s %18s %18s %16s\n', ...
        'p','CoreStorage','RelativeError','TotalStorage');

    for t = 1:count
        p = pValues(t);
        [model,qrStats] = rc_tlr_build(CP,grids,opts.domain,p,opts.qrTol,opts.strength);
        ranks(t,:) = model.ranks;
        coreStorage{t} = sprintf('%d x %d x %d',model.ranks);
        C = model.C;

        % The factored representation is complete before this evaluation.
        X = model.U*model.E1;
        Y = model.V*model.E2;
        Z = model.W*model.E3;
        totalStorage(t) = numel(C)+numel(X)+numel(Y)+numel(Z);
        [relativeError(t),metrics] = rc_htlr_relative_kernel_error( ...
            C,X,Y,Z,CP,opts.validation,opts.blockSize);
        diagnostics{t} = struct('qr',{qrStats},'metrics',metrics);
        fprintf('%5d %18s %18.8e %16.0f\n', ...
            p,coreStorage{t},relativeError(t),totalStorage(t));
        if opts.stopWhenMet && relativeError(t)<=opts.tol, break; end
    end
    take = 1:t;
    diagnostics = diagnostics(take);
    results = table(pValues(take),ranks(take,1),ranks(take,2),ranks(take,3), ...
        coreStorage(take),relativeError(take), ...
        totalStorage(take),relativeError(take)<=opts.tol, ...
        'VariableNames',{'p','rX','rY','rZ','CoreStorage','RelativeError', ...
        'TotalStorage','ToleranceMet'});
end
