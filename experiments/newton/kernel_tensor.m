function [K, x, y, z, RS] = kernel_tensor(n, outputFormat, part, dataDir)
%KERNEL_TENSOR Authors' range-separated Newton CP benchmark on [-1,1]^3.
%   [K,x,y,z] = kernel_tensor(n) returns the dense LONG-RANGE tensor.
%   [CP,x,y,z,RS] = kernel_tensor(n,'cp') avoids forming the n^3 tensor.
%   [S,x,y,z,RS] = kernel_tensor(n,'slice') returns the middle z-slice.
%   Optional PART: 'long' (default), 'short', or 'full'.
%   Optional DATADIR: directory containing Newt_canon_<n>.mat.
%
%   CP.lambda and CP.U{1:3} define sum_s lambda(s)*u_s*v_s*w_s.
%   RS.long and RS.short contain both parts; RS.long + RS.short = full CP.
%   RS.selectedColumns identifies the original data columns used for K.
%   In 'slice' mode z is the scalar slice coordinate; RS.grid is the full grid.
%
%   Matches bonans/ChebTuck/newton_tabs.m, commit
%   2aa3729fe9991377bf1a6815ab3edd73423689dd:
%       h = 2/n; U = U1r/h; Rl = ceil(R/2)-3;
%       long columns = 1:Rl; grid = cell centres.
%   The supplied factors approximate the cell-AVERAGED Newton kernel.
%   They are not point samples of 1/sqrt(x^2+y^2+z^2).
%   No artificial regularization, new quadrature, or CP recompression is used.
%   The dataset covers n = 256,512,1024,2048,4096,8192.
%   Dense output uses 8*n^3 bytes: 4 TiB at n=8192.
%   Use 'cp' or 'slice' for large n.
%   No Tensor Toolbox or Chebfun dependency.

    validateattributes(n, {'numeric'}, ...
        {'scalar','real','finite','integer','positive'}, mfilename, 'n');
    n = double(n);
    if ~ismember(n, [256 512 1024 2048 4096 8192])
        error('kernel_tensor:UnsupportedGrid', ...
            'Original data support n = 256,512,1024,2048,4096,8192.');
    end
    if nargin < 2 || isempty(outputFormat), outputFormat = 'dense'; end
    if nargin < 3 || isempty(part), part = 'long'; end
    outputFormat = validatestring(outputFormat, {'dense','cp','slice'});
    part = validatestring(part, {'long','short','full'});
    if nargin < 4 || isempty(dataDir)
        dataDir = fullfile(fileparts(mfilename('fullpath')), ...
            'data', 'Newton_CP_-1_1');
    end
    dataFile = fullfile(dataDir, sprintf('Newt_canon_%d.mat', n));
    if ~isfile(dataFile)
        error('kernel_tensor:MissingData', ...
            'Missing %s. Extract the complete supplied ZIP, including data.', dataFile);
    end
    S = load(dataFile, 'U1r');
    if ~isfield(S, 'U1r')
        error('kernel_tensor:InvalidData', 'The data file must contain U1r.');
    end
    validateattributes(S.U1r, {'double'}, ...
        {'2d','real','finite','nonempty','nrows',n}, mfilename, 'U1r');
    h = 2/n;
    U = S.U1r / h;
    R = size(U, 2);
    Rl = ceil(R/2)-3;
    if Rl < 1 || Rl >= R
        error('kernel_tensor:InvalidRank', 'Data rank is incompatible with the benchmark split.');
    end
    x = linspace(-1+h/2, 1-h/2, n);
    y = x;
    z = x;
    switch part
        case 'long', indices = 1:Rl;
        case 'short', indices = Rl+1:R;
        case 'full', indices = 1:R;
    end
    V = U(:, indices);
    CP = struct('lambda', ones(numel(indices),1), 'U', {{V,V,V}});
    switch outputFormat
        case 'cp'
            K = CP;
        case 'slice'
            k0 = floor(n/2);
            K = (V .* V(k0,:)) * V.';
            z = x(k0);
        case 'dense'
            K = zeros(n,n,n);
            Vt = V.';
            for k = 1:n
                K(:,:,k) = (V .* V(k,:)) * Vt;
            end
    end
    if nargout >= 5
        L = U(:,1:Rl);
        T = U(:,Rl+1:R);
        RS = struct;
        RS.long = struct('lambda',ones(Rl,1),'U',{{L,L,L}});
        RS.short = struct('lambda',ones(R-Rl,1),'U',{{T,T,T}});
        RS.fullRank = R;
        RS.longRank = Rl;
        RS.shortRank = R-Rl;
        RS.h = h;
        RS.grid = x;
        RS.domain = repmat([-1 1],3,1);
        RS.cellEdges = linspace(-1,1,n+1);
        RS.middleIndex = floor(n/2);
        RS.selectedPart = part;
        RS.selectedColumns = indices;
        RS.sourceFile = dataFile;
        RS.sourceCommit = '2aa3729fe9991377bf1a6815ab3edd73423689dd';
    end
end
