function C = selected_core(p,I,F)
%SELECTED_CORE Extract unweighted skeleton from sampled CP factors or dense G.
%   C=selected_core(p,I,G), where G is output 5 of weighted_cheb_tensor.
    assert(nargin>=3 && ~isempty(F),'Supply unweighted samples.');
    if isscalar(p), p = repmat(p,1,3); end
    assert(iscell(I) && numel(I)==3,'I must contain three index lists.');
    for d = 1:3
        validateattributes(I{d},{'numeric'}, ...
            {'vector','integer','positive','<=',p(d),'nonempty'});
        assert(numel(unique(I{d}))==numel(I{d}),'Repeated selected index.');
    end
    if isstruct(F)
        assert(~isfield(F,'weighted') || ~F.weighted, ...
            'Extract the core from unweighted samples, not Gw.');
        assert(isequal(cellfun(@(u) size(u,1),F.U),reshape(p,1,3)), ...
            'p does not match the sampled factors.');
        C = rc_tlr_cp_values(F,I);
    else
        C = F(I{1},I{2},I{3});
    end
end
