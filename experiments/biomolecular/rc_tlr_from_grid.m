function [M,G] = rc_tlr_from_grid(F,grids,p,tol,options)
%RC_TLR_FROM_GRID Algorithm 2: full-grid cubic spline -> RC functional Alg. 1.
% [M,G]=rc_tlr_from_grid(F,{x,y,z},129,1e-9);
% F uses ndgrid ordering, size n1-by-n2-by-n3; each grid has >=4 points.
% Spline construction and sampling count in TotalTime. Input preparation,
% evaluated output factors, validation and tensor reconstruction are excluded.
    if nargin<5, options=struct; end
    validateattributes(F,{'double','single'},{'real','finite','nonempty'});
    assert(ndims(F)==3 && iscell(grids) && numel(grids)==3,'Supply a full 3-D tensor and its three grids.');
    for d=1:3
        g=grids{d}(:);
        assert(numel(g)==size(F,d) && numel(g)>=4 && all(isfinite(g)) && all(diff(g)>0), ...
            'Each grid must match F and have at least four increasing points.');
        assert(g(1)>=-1 && g(end)<=1,'Input grids must lie within [-1,1].');
        grids{d}=g;
    end
    t=tic;
    Q=griddedInterpolant(grids,double(F),'spline','spline');
    [M,G]=rc_tlr_functional(@(x,y,z)Q(x,y,z),p,tol,options);
    M.InputSize=[size(F,1),size(F,2),size(F,3)];
    M.InputRepresentation='full_grid_spline';
    M.BuildSeconds=toc(t); M.TotalTime=M.BuildSeconds;
end
