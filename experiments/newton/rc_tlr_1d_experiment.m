function [results, model, diagnostics] = rc_tlr_1d_experiment(options)
%RC_TLR_1D_EXPERIMENT RC-TLR on the long-range Newton CP input.
%   The defaults are n=256, pValues=[129 257 513 1025], relative-error
%   target tol=1e-7, and weighted RRQR tolerance qrTol=5e-8. The default
%   validation metric is relative maximum error on the whole middle slice.
%   Set validation='full' for the full-volume relative Frobenius error.
%
%   BuildSeconds is the median of repeated constructions of C,E1,E2,E3,U,V,W.
%   TotalTime is the same value for callers using the previous column name.
%   Data loading, U*E1, V*E2, W*E3, reconstruction and validation are outside
%   this timer. TotalStorage counts the core and three evaluated grid factors.
%
%   Defaults: warmups=1, repeats=5, strength=2, blockSize=8, stopWhenMet=true.
%   For custom input, provide CP, grids={x,y,z}, and the physical domain.
    if nargin < 1 || isempty(options), options = struct; end
    assert(isstruct(options) && isscalar(options),'Options must be a scalar struct.');
    opts = struct('n',256,'pValues',[129 257 513 1025], ...
        'tol',1e-7,'qrTol',5e-8,'strength',2,'blockSize',8, ...
        'warmups',1,'repeats',5,'validation','slice','stopWhenMet',true, ...
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
    validateattributes(opts.warmups,{'numeric'},{'scalar','integer','nonnegative','finite'});
    validateattributes(opts.repeats,{'numeric'},{'scalar','integer','positive','finite'});
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
    buildSeconds = zeros(count,1);
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
    fprintf('BuildSeconds includes construction of the selected core and six matrices.\n');
    fprintf('%5s %18s %18s %18s %16s\n', ...
        'p','CoreStorage','RelativeError','BuildSeconds (s)','TotalStorage');

    for t = 1:count
        p = pValues(t);
        times = zeros(opts.repeats,1);
        for rep = 1:opts.warmups+opts.repeats
            timer = tic;
            [model,qrStats] = rc_tlr_build(CP,grids,opts.domain,p,opts.qrTol,opts.strength);
            seconds = toc(timer);
            if rep > opts.warmups
                times(rep-opts.warmups) = seconds;
            end
        end
        buildSeconds(t) = median(times);
        ranks(t,:) = model.ranks;
        coreStorage{t} = sprintf('%d x %d x %d',model.ranks);
        C = model.C;

        % The factored representation is complete before this evaluation.
        X = model.U*model.E1;
        Y = model.V*model.E2;
        Z = model.W*model.E3;
        totalStorage(t) = numel(C)+numel(X)+numel(Y)+numel(Z);
        [relativeError(t),~,metrics] = rc_htlr_relative_kernel_error( ...
            C,X,Y,Z,CP,opts.validation,opts.blockSize);
        diagnostics{t} = struct('qr',{qrStats},'metrics',metrics,'buildRepeats',times);
        fprintf('%5d %18s %18.8e %18.6f %16.0f\n', ...
            p,coreStorage{t},relativeError(t),buildSeconds(t),totalStorage(t));
        if opts.stopWhenMet && relativeError(t)<=opts.tol, break; end
    end
    take = 1:t;
    diagnostics = diagnostics(take);
    results = table(pValues(take),ranks(take,1),ranks(take,2),ranks(take,3), ...
        coreStorage(take),relativeError(take),buildSeconds(take),buildSeconds(take), ...
        totalStorage(take),relativeError(take)<=opts.tol, ...
        'VariableNames',{'p','rX','rY','rZ','CoreStorage','RelativeError', ...
        'BuildSeconds','TotalTime','TotalStorage','ToleranceMet'});
end
