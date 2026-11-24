function F = rc_htlr_cp_values(CP, indices)
%RC_HTLR_CP_VALUES Evaluate a Cartesian block of a real CP tensor.
%   CP.lambda is R-by-1; CP.U is a cell of three n_d-by-R factors.
    if nargin < 2
        indices = {1:size(CP.U{1},1),1:size(CP.U{2},1),1:size(CP.U{3},1)};
    end
    A = CP.U{1}(indices{1},:);
    B = CP.U{2}(indices{2},:);
    D = CP.U{3}(indices{3},:);
    A = A .* CP.lambda(:).';
    Bt = B.';
    F = zeros(size(A,1),size(B,1),size(D,1));
    for k = 1:size(D,1)
        F(:,:,k) = (A .* D(k,:))*Bt;
    end
end
