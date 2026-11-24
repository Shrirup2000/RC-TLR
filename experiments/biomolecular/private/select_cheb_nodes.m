function [I,selectedNodes,ranks,info] = select_cheb_nodes(F_w,tol,strength)
%SELECT_CHEB_NODES Independent strong RRQR in each weighted tensor mode.
    if nargin < 2 || isempty(tol), tol = 1e-9; end
    if nargin < 3 || isempty(strength), strength = 2; end
    I = cell(1,3); selectedNodes = cell(1,3); info = cell(1,3);
    ranks = zeros(1,3);
    for d = 1:3
        [I{d},selectedNodes{d},ranks(d),info{d}] = ...
            strong_rrqr_nodes(F_w,d,tol,strength);
    end
end
