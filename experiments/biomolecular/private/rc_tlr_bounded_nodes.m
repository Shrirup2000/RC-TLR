function [I,r,info] = rc_tlr_bounded_nodes(Gw,tol,strength,workspaceMB)
%RC_TLR_BOUNDED_NODES Exact CP reduction followed by the supplied strong RRQR.
% Same weighted unfolding as strong_rrqr_nodes; its large Khatri-Rao array
% is reduced by row blocks. No Gram matrix, SVD, or tensor truncation is used.
% workspaceMB controls block sizing, not measured peak process RAM.
    if nargin<4, workspaceMB=128; end
    I=cell(1,3); info=cell(1,3); r=zeros(1,3);
    R=numel(Gw.lambda);
    for d=1:3
        other=setdiff(1:3,d,'stable');
        [~,B]=qr(Gw.U{other(1)},0);
        [~,D]=qr(Gw.U{other(2)},0);
        block=max(1,floor(workspaceMB*2^20/(4*8*R)));
        S=zeros(0,R); nr=size(B,1)*size(D,1);
        for first=1:block:nr
            ids=first:min(first+block-1,nr);
            [i,j]=ind2sub([size(B,1),size(D,1)],ids);
            P=B(i,:).*D(j,:);
            [~,S]=qr([S;P],0);
        end
        A=(S.*Gw.lambda(:).')*Gw.U{d}.';
        p=size(Gw.U{d},1);
        % Feed the exact small unfolding to the supplied RRQR implementation.
        [I{d},~,r(d),info{d}]=strong_rrqr_nodes( ...
            reshape(A.',p,size(A,1),1),1,tol,strength);
        info{d}.mode=d;
        info{d}.weights=Gw.weights{d};
        info{d}.route='CP block orthogonal reduction';
        assert(info{d}.accuracySatisfied && info{d}.strongSatisfied, ...
            'RC-TLR QR accuracy/strength check failed in mode %d.',d);
    end
end
