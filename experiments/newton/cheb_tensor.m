function [F,x_cheb,y_cheb,z_cheb] = cheb_tensor(p,CP,grids,domain,representation)
%CHEB_TENSOR Spline-sample CP factors at first-kind Chebyshev nodes.
%   F = cheb_tensor(p,CP,{x,y,z},domain) returns CP form, not a dense array.
%   Set representation='dense' only for verification on small problems.
%   No arguments after p: use the n=256 long-range Newton input.
%   Splines extrapolate from cell centres to the ends of the physical box.
    if nargin < 2 || isempty(CP)
        [CP,x,y,z] = kernel_tensor(256,'cp'); grids = {x,y,z};
    elseif nargin < 3 || isempty(grids)
        error('cheb_tensor:MissingGrid','Provide the original CP factor grids.');
    end
    if nargin < 4 || isempty(domain), domain = repmat([-1 1],3,1); end
    if nargin < 5 || isempty(representation), representation = 'cp'; end
    representation = validatestring(representation,{'cp','dense'});
    validateattributes(p,{'numeric'},{'vector','finite','integer','positive'});
    if isscalar(p), p = repmat(p,1,3); end
    assert(numel(p)==3,'p must have one or three entries.');
    validateattributes(domain,{'numeric'},{'size',[3 2],'real','finite'});
    assert(all(domain(:,2)>domain(:,1)),'Invalid domain endpoints.');
    assert(isstruct(CP) && isfield(CP,'U') && isfield(CP,'lambda'), ...
        'CP must contain U and lambda.');
    assert(iscell(CP.U) && numel(CP.U)==3 && iscell(grids) && numel(grids)==3, ...
        'Supply three CP factors and three coordinate grids.');
    validateattributes(CP.lambda,{'numeric'},{'vector','real','finite','nonempty'});
    R = numel(CP.lambda); P = cell(1,3); nodes = cell(1,3);
    for d = 1:3
        validateattributes(CP.U{d},{'numeric'},{'2d','real','finite','ncols',R});
        g = double(grids{d}(:));
        assert(numel(g)==size(CP.U{d},1) && numel(g)>=4 && ...
            all(isfinite(g)) && all(diff(g)>0), ...
            'Each CP grid must have at least four strictly increasing coordinates.');
        assert(g(1)>=domain(d,1) && g(end)<=domain(d,2), ...
            'The source grid must lie within its physical interval.');
        nodes{d} = domain(d,1)+(domain(d,2)-domain(d,1))*(my_chev2(p(d))+1)/2;
        P{d} = interp1(g,double(CP.U{d}),nodes{d}(:),'spline','extrap');
    end
    F = struct('lambda',double(CP.lambda(:)),'U',{P}, ...
        'nodes',{nodes},'domain',domain,'weighted',false);
    x_cheb = nodes{1}; y_cheb = nodes{2}; z_cheb = nodes{3};
    if strcmp(representation,'dense'), F = rc_htlr_cp_values(F); end
end
