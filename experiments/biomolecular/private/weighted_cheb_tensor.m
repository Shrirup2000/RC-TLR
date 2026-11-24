function [F_w,x_cheb,y_cheb,z_cheb,F] = weighted_cheb_tensor(p,varargin)
%WEIGHTED_CHEB_TENSOR Sample once, return weighted and original CP tensors.
%   Arguments after p are forwarded to cheb_tensor.
    [F,x_cheb,y_cheb,z_cheb] = cheb_tensor(p,varargin{:});
    [WX,WY,WZ] = cheb_weight_matrices(p);
    F_w = weighted_core(F,WX,WY,WZ);
end
