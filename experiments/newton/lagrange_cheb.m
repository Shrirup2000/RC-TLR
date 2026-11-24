function L = lagrange_cheb(xi,x,baryWeights)
%LAGRANGE_CHEB Barycentric Lagrange evaluation (full p-node basis in RC-HTLR).
%   Optional exact barycentric weights avoid quadratic weight construction.
    xi = xi(:).'; x = x(:); p = numel(xi);
    assert(p>0 && all(isfinite(xi)) && all(isfinite(x)), 'Invalid nodes or points.');
    assert(numel(unique(xi))==p,'Interpolation nodes must be distinct.');
    if p==1, L = ones(numel(x),1); return; end
    if nargin < 3 || isempty(baryWeights)
        % Fast path for the full first-kind reference grid.
        if isequal(xi,my_chev2(p))
            baryWeights = (-1).^(0:p-1).*sin((2*(1:p)-1)*pi/(2*p));
        else
            D = xi.'-xi; D(1:p+1:end) = 1;
            logs = -sum(log(abs(D)),2);
            baryWeights = (prod(sign(D),2).*exp(logs-max(logs))).';
        end
    end
    w = baryWeights(:).';
    assert(numel(w)==p && all(isfinite(w)) && all(w~=0),'Invalid barycentric weights.');
    w = w/max(abs(w));
    L = zeros(numel(x),p);
    [hit,index] = ismember(x,xi); rows = find(hit);
    L(sub2ind(size(L),rows,index(rows))) = 1;
    rows = find(~hit);
    if ~isempty(rows)
        D = x(rows)-xi;
        T = (min(abs(D),[],2)./D).*w;
        denom = sum(T,2);
        assert(all(isfinite(denom)) && all(denom~=0),'Unstable barycentric evaluation.');
        L(rows,:) = T./denom;
    end
end
