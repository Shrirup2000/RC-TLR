function F_w = weighted_core(F,WX,WY,WZ)
%WEIGHTED_CORE Square-root weighting of a dense tensor or its CP factors.
    weights = {WX(:),WY(:),WZ(:)};
    for d = 1:3
        validateattributes(weights{d},{'numeric'},{'real','finite','positive'});
    end
    if isstruct(F)
        if isfield(F,'weighted') && F.weighted
            error('weighted_core:AlreadyWeighted','Do not apply weights twice.');
        end
        F_w = F;
        for d = 1:3
            assert(numel(weights{d})==size(F.U{d},1),'Weight length mismatch.');
            F_w.U{d} = F.U{d} .* sqrt(weights{d});
        end
        F_w.weights = weights;
        F_w.weighted = true;
    else
        assert(isequal([size(F,1),size(F,2),size(F,3)],cellfun(@numel,weights)), ...
            'Weight lengths must match tensor dimensions.');
        F_w = F .* reshape(sqrt(WX),[],1,1) .* ...
            reshape(sqrt(WY),1,[],1) .* reshape(sqrt(WZ),1,1,[]);
    end
end
