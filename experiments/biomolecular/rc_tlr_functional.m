function [M,G] = rc_tlr_functional(fun,p,tol,options)
%RC_TLR_FUNCTIONAL Algorithm 1: weighted strong RRQR of FUNCTION VALUES.
% [M,G]=rc_tlr_functional(@(x,y,z)exp(-(x.*y.*z).^2),[17 19 21],1e-9);
% M.C is G(Ix,Iy,Iz); M.H contains unweighted p_d-by-r_d interpolation maps.
% M.TotalTime=M.BuildSeconds: construction of the continuous approximation.
% rc_tlr_functional_factors(M,{x,y,z}) evaluates Lagrange basis times M.H.
% G is a diagnostic output; M does not retain G, the spline, or the input.
    if nargin<4, options=struct; end
    c=struct('strength',2,'maxSampleMB',128);
    names=fieldnames(options);
    for i=1:numel(names)
        assert(isfield(c,names{i}),'Unknown functional option: %s',names{i});
        c.(names{i})=options.(names{i});
    end
    validateattributes(p,{'numeric'},{'vector','integer','positive','finite'});
    if isscalar(p), p=repmat(p,1,3); end
    assert(numel(p)==3 && isa(fun,'function_handle'),'Supply a function handle and three node counts.');
    p=reshape(double(p),1,3);
    validateattributes(tol,{'numeric'},{'scalar','real','>',0,'<',1,'finite'});
    validateattributes(c.strength,{'numeric'},{'scalar','>',1,'finite'});
    validateattributes(c.maxSampleMB,{'numeric'},{'scalar','positive','finite'});
    assert(8*prod(p)<=c.maxSampleMB*2^20, ...
        'Dense G exceeds maxSampleMB. Reduce p or explicitly raise the functional sampling limit.');
    t=tic;
    domain=repmat([-1 1],3,1);
    [Gw,x,y,z,G]=weighted_cheb_tensor(p,fun,[],domain,'dense');
    [I,~,r,q]=select_cheb_nodes(Gw,tol,c.strength);
    for d=1:3
        assert(q{d}.accuracySatisfied && q{d}.strongSatisfied, ...
            'Functional RC-TLR QR check failed in mode %d.',d);
    end
    C=selected_core(p,I,G);
    [~,~,~,H]=selected_factor_matrices({[],[],[]},p,I,q,domain);
    M=struct('C',C,'H',{H},'nodes',{{x,y,z}},'p',p,'Indices',{I}, ...
        'Ranks',r,'Domain',domain,'Tolerance',tol,'Strength',c.strength, ...
        'MaxQRResidual',max(cellfun(@(s)s.relativeResidual,q)), ...
        'MaxRho',max(cellfun(@(s)s.maxRho,q)), ...
        'ContinuousStorage',numel(C)+sum(cellfun(@numel,H)), ...
        'SampleEntries',prod(p),'BuildSeconds',NaN,'TotalTime',NaN);
    M.BuildSeconds=toc(t); M.TotalTime=M.BuildSeconds;
end
