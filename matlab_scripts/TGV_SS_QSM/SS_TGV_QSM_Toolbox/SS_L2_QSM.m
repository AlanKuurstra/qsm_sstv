function out = SS_L2_QSM(params)

    pcg_tol = 1e-3;
    pcg_maxit = 100;
    alpha = params.alpha;
    N = params.N;
    M = params.M;
    H = params.H;
    D = params.D;
    numKernels = size(H,4);

    % For gradient
    [k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);
    E1 = 1 - exp(2i .* pi .* k1 / N(1));
    E2 = 1 - exp(2i .* pi .* k2 / N(2));
    E3 = 1 - exp(2i .* pi .* k3 / N(3));
    EtE = abs(E1).^2 + abs(E2).^2 + abs(E3).^2;

    % Precompute
    F_phi = fftn(params.phase_unwrap);
    HD = bsxfun(@times, H, D);
    DtHt = conj(HD);
    MtM = abs(M).^2;
    
    % Compute rhs
    rhs = zeros(N);
    for h = 1:numKernels
        rhs = rhs + DtHt(:,:,:,h) .* fftn(MtM(:,:,:,h) .* ifftn(H(:,:,:,h) .* F_phi));
    end

    % Precondition
    DtHtMtMHD = sum(DtHt.*MtM.*HD,4);
    precond = 1 ./ (eps + DtHtMtMHD(:) + alpha * EtE(:));
    pre_inv = @(x, precond) precond .* x;
    
    tic
    Fx = pcg(@(Fx)apply_SS_L2( Fx, HD, DtHt, MtM, EtE, alpha, N, numKernels), rhs(:), pcg_tol, pcg_maxit, @(x)pre_inv(x, precond));
    toc
    
    out.x = real(ifftn(reshape(Fx, N)));

end

function out = apply_SS_L2( Fx, HD, DtHt, MtM, EtE, alpha, N, numKernels)

    Fx = reshape(Fx, N);
    res = alpha * EtE .*  Fx;

    for h = 1:numKernels
        res = res + DtHt(:,:,:,h).*fftn(MtM(:,:,:,h).* ifftn( HD(:,:,:,h).*Fx));   
    end

    out = res(:);

end