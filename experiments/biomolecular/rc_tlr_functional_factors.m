function [A,L] = rc_tlr_functional_factors(M,grids)
%RC_TLR_FUNCTIONAL_FACTORS A_d=L_d*E_d, for Algorithm 1's final grid factors.
% [A,L]=rc_tlr_functional_factors(M,{x,y,z}) also returns the uncompressed
% Lagrange matrices L{d}, of size n_d-by-p_d. A{d} is n_d-by-r_d.
% Evaluation is outside the continuous-model construction timer.
    assert(iscell(grids) && numel(grids)==3,'Supply three evaluation vectors.');
    A=cell(1,3); L=cell(1,3);
    for d=1:3
        x=grids{d}(:); validateattributes(x,{'numeric'},{'real','finite'});
        A{d}=zeros(numel(x),size(M.H{d},2));
        if nargout>1, L{d}=zeros(numel(x),M.p(d)); end
        for first=1:128:numel(x)
            rows=first:min(first+127,numel(x));
            basis=lagrange_cheb(M.nodes{d},x(rows));
            A{d}(rows,:)=basis*M.H{d};
            if nargout>1, L{d}(rows,:)=basis; end
        end
    end
end
