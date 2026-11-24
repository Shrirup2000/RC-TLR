function [WX,WY,WZ,DX,DY,DZ] = cheb_weight_matrices(p)
%CHEB_WEIGHT_MATRICES Positive reference weights, scalar p or [px py pz].
    validateattributes(p, {'numeric'}, {'vector','integer','positive','finite'});
    if isscalar(p), p = repmat(p,1,3); end
    assert(numel(p)==3, 'p must be scalar or contain three node counts.');
    weights = cell(1,3);
    for d = 1:3
        k = (1:p(d)).';
        weights{d} = (pi/p(d))*sin((2*k-1)*pi/(2*p(d)));
    end
    WX = weights{1}; WY = weights{2}; WZ = weights{3};
    if nargout > 3
        DX = spdiags(WX,0,p(1),p(1));
        DY = spdiags(WY,0,p(2),p(2));
        DZ = spdiags(WZ,0,p(3),p(3));
    end
end
