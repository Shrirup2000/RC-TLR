function [U,V,W,H] = selected_factor_matrices(grids,p,I,info,domain)
%SELECTED_FACTOR_MATRICES Undo weights, then evaluate the full Lagrange basis.
%   GRIDS={x,y,z} uses actual physical coordinates (preferred).
%   A scalar/vector of grid sizes instead creates cell-centred grids on DOMAIN.
%   Optional fourth output H contains the continuous p_d-by-r_d maps.
    if nargin < 5 || isempty(domain), domain = repmat([-1 1],3,1); end
    if isscalar(p), p = repmat(p,1,3); end
    if ~iscell(grids)
        n = grids; if isscalar(n), n = repmat(n,1,3); end
        validateattributes(n,{'numeric'},{'numel',3,'integer','positive','finite'});
        grids = cell(1,3);
        for d = 1:3
            h = (domain(d,2)-domain(d,1))/n(d);
            grids{d} = linspace(domain(d,1)+h/2,domain(d,2)-h/2,n(d));
        end
    end
    factors = cell(1,3); H = cell(1,3);
    for d = 1:3
        pd = p(d); data = info{d};
        assert(data.p==pd && isequal(I{d}(:),reshape(data.piv(1:data.rank),[],1)), ...
            'Selected indices do not match their QR data.');
        if isempty(data.weights)
            theta = (2*(1:pd).'-1)*pi/(2*pd);
            weights = (pi/pd)*sin(theta);
        else
            weights = data.weights(:);
        end
        s = sqrt(weights);
        H{d} = data.E.' .* (s(I{d}).'./s);
        % Exact cardinal rows also avoid small roundoff in the scaling ratio.
        H{d}(I{d},:) = eye(data.rank);
        xref = 2*(grids{d}(:)-domain(d,1))/(domain(d,2)-domain(d,1))-1;
        nodes = my_chev2(pd);
        bary = (-1).^(0:pd-1).*sin((2*(1:pd)-1)*pi/(2*pd));
        factors{d} = zeros(numel(xref),data.rank);
        for first = 1:256:numel(xref)
            rows = first:min(first+255,numel(xref));
            B = lagrange_cheb(nodes,xref(rows),bary);
            factors{d}(rows,:) = B*H{d};
        end
    end
    U = factors{1}; V = factors{2}; W = factors{3};
end
