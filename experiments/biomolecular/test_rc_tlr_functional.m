function test_rc_tlr_functional()
%TEST_RC_TLR_FUNCTIONAL Algorithms 1/2; requires only base MATLAB.
    f=@(x,y,z) .7+2*x-3*y+.2*(2*x.^2-1).*y.*(4*z.^3-3*z);
    p=[9 11 13]; tol=1e-11;
    [M,G]=rc_tlr_functional(f,p,tol);
    g={linspace(-1,1,15),linspace(-.9,.8,17),linspace(-1,1,19)};
    [A,L]=rc_tlr_functional_factors(M,g);
    for d=1:3
        assert(isequal(size(L{d}),[numel(g{d}),p(d)]));
        assert(norm(A{d}-L{d}*M.H{d},'fro')<1e-12*max(1,norm(A{d},'fro')));
        assert(isequal(M.H{d}(M.Indices{d},:),eye(M.Ranks(d))));
    end
    k=6; Y=slice(M.C,A,k); [X,Q]=ndgrid(g{1},g{2});
    ref=f(X,Q,g{3}(k));
    assert(norm(Y-ref,'fro')<1e-9*norm(ref,'fro'));
    check_samples(M,G,f); assert(M.TotalTime==M.BuildSeconds);
    assert(M.ContinuousStorage==prod(M.Ranks)+sum(p.*M.Ranks));
    assert(~isfield(M,'G') && ~isfield(M,'Q') && ~isfield(M,'F'));
    [B,~]=rc_tlr_functional(@(x,y,z)0*x+2,1,tol);
    assert(B.C==2 && all(B.Ranks==1));
    [B,~]=rc_tlr_functional(@(x,y,z)0*x,[3 4 5],tol);
    assert(all(B.C(:)==0));
    % Non-cubic, asymmetric full tensor: no CP factors enter Algorithm 2.
    n=[13 14 15]; source=arrayfun(@(a)linspace(-1+1/a,1-1/a,a),n,'UniformOutput',false);
    [X,Y,Z]=ndgrid(source{:}); F=f(X,Y,Z);
    [B,G]=rc_tlr_from_grid(F,source,[7 9 11],tol);
    check_samples(B,G,f); assert(B.TotalTime==B.BuildSeconds);
    assert(isequal(B.InputSize,n));
    [A,~]=rc_tlr_functional_factors(B,source);
    Y=slice(B.C,A,7); ref=F(:,:,7);
    assert(norm(Y-ref,'fro')<1e-9*norm(ref,'fro'));
    % General nonseparable data: independently compare spline samples.
    F=sin(X.*Y+2*Z)+.1*cos(3*X+Y.*Z);
    [B,G]=rc_tlr_from_grid(F,source,[7 9 11],tol);
    [X,Y,Z]=ndgrid(B.nodes{:});
    ref=interpn(source{1},source{2},source{3},F,X,Y,Z,'spline');
    assert(norm(G(:)-ref(:))<1e-11*norm(ref(:)));
    fprintf('RC functional and full-grid spline tests passed.\n');
end

function check_samples(M,G,f)
    [X,Y,Z]=ndgrid(M.nodes{:}); ref=f(X,Y,Z);
    assert(norm(G(:)-ref(:))<1e-11*max(norm(ref(:)),1));
    assert(isequal(M.C,G(M.Indices{1},M.Indices{2},M.Indices{3})));
    [A,~]=rc_tlr_functional_factors(M,M.nodes);
    for k=1:M.p(3)
        Y=slice(M.C,A,k); ref=G(:,:,k);
        assert(norm(Y-ref,'fro')<1e-9*max(norm(ref,'fro'),1));
    end
end

function Y=slice(C,A,k)
    r=[size(C,1),size(C,2),size(C,3)];
    small=reshape(reshape(C,r(1)*r(2),r(3))*A{3}(k,:).',r(1),r(2));
    Y=(A{1}*small)*A{2}.';
end
