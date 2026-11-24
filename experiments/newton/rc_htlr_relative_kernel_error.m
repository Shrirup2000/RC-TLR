function [relativeError,metrics] = rc_htlr_relative_kernel_error( ...
    C,U,V,W,CP,validation,blockSize)
%RC_HTLR_RELATIVE_KERNEL_ERROR Direct CP reference; no full reference tensor.
%   'slice': paper's relative maximum error on k=floor(nz/2).
%   'full' : relative Frobenius error over all entries.
    if nargin < 6 || isempty(validation), validation = 'slice'; end
    if nargin < 7 || isempty(blockSize), blockSize = 8; end
    validation = validatestring(validation,{'slice','full'});
    validateattributes(blockSize,{'numeric'},{'scalar','integer','positive','finite'});
    nx = size(U,1); ny = size(V,1); nz = size(W,1);
    rx = size(U,2); ry = size(V,2); rz = size(W,2);
    assert(isequal(cellfun(@(a) size(a,1),CP.U),[nx ny nz]),'Reference grid mismatch.');
    assert(isequal([size(C,1) size(C,2) size(C,3)],[rx ry rz]),'Core size mismatch.');
    errorSquared = 0; referenceSquared = 0; errorMax = 0; referenceMax = 0;
    if strcmp(validation,'slice')
        k0 = max(1,floor(nz/2));
        sliceCore = reshape(reshape(C,rx*ry,rz)*W(k0,:).',rx,ry);
        for first = 1:blockSize:nx
            ix = first:min(first+blockSize-1,nx);
            approximation = (U(ix,:)*sliceCore)*V.';
            reference = rc_htlr_cp_values(CP,{ix,1:ny,k0});
            delta = approximation-reference;
            errorSquared = errorSquared+sum(abs(delta(:)).^2);
            referenceSquared = referenceSquared+sum(abs(reference(:)).^2);
            errorMax = max(errorMax,max(abs(delta(:))));
            referenceMax = max(referenceMax,max(abs(reference(:))));
        end
    else
        % Tile x as well as z: avoid a persistent nx-by-ny-by-rz intermediate.
        for firstX = 1:blockSize:nx
            ix = firstX:min(firstX+blockSize-1,nx); bx = numel(ix);
            H = reshape(U(ix,:)*reshape(C,rx,ry*rz),bx,ry,rz);
            H = V*reshape(permute(H,[2 1 3]),ry,bx*rz);
            H = permute(reshape(H,ny,bx,rz),[2 1 3]);
            H = reshape(H,bx*ny,rz);
            for firstZ = 1:blockSize:nz
                iz = firstZ:min(firstZ+blockSize-1,nz);
                approximation = reshape(H*W(iz,:).',bx,ny,numel(iz));
                reference = rc_htlr_cp_values(CP,{ix,1:ny,iz});
                delta = approximation-reference;
                errorSquared = errorSquared+sum(abs(delta(:)).^2);
                referenceSquared = referenceSquared+sum(abs(reference(:)).^2);
                errorMax = max(errorMax,max(abs(delta(:))));
                referenceMax = max(referenceMax,max(abs(reference(:))));
            end
        end
    end
    frob = safe_ratio(sqrt(errorSquared),sqrt(referenceSquared));
    maxerr = safe_ratio(errorMax,referenceMax);
    if strcmp(validation,'slice'), relativeError = maxerr;
    else, relativeError = frob; end
    metrics = struct('scope',validation,'relativeFrobenius',frob, ...
        'relativeMaximum',maxerr);
end

function value = safe_ratio(numerator,denominator)
    if denominator==0
        if numerator==0, value=0; else, value=Inf; end
    else
        value = numerator/denominator;
    end
end
