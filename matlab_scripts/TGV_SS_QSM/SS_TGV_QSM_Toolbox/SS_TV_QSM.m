
function out = SS_TV_QSM(params)

    % Retrieve data
    alpha = params.alpha;
    mu1 = params.mu1;
    mu2 = params.mu2;
    maxOuterIter = params.maxOuterIter;
    tol_soln = params.tol_soln;
    N = params.N;
    M = params.M;
    H = params.H;
    D = params.D;
    numKernels = size(H,4);

    %% Precompute
    
    % For objective function
    F_phi = fftn(params.phase_unwrap);
    MFtHF_phi = zeros([N,numKernels]);
    for i = 1:numKernels
        MFtHF_phi(:,:,:,i) = M(:,:,:,i).*ifftn(H(:,:,:,i).*F_phi);
    end
    
    % For others
    [k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);
    E1 = 1 - exp(2i .* pi .* k1 / N(1));
    E2 = 1 - exp(2i .* pi .* k2 / N(2));
    E3 = 1 - exp(2i .* pi .* k3 / N(3));

    E1t = conj(E1); E2t = conj(E2); E3t = conj(E3);
    EtE = abs(E1).^2 + abs(E2).^2 + abs(E3).^2;
 
    Dt = conj(D); Ht = conj(H); Mt = conj(M);
    MtM = abs(M).^2;
    MtM_FtHF_phi = zeros([N,numKernels]);
    for i = 1:numKernels
        MtM_FtHF_phi(:,:,:,i) = Mt(:,:,:,i).*MFtHF_phi(:,:,:,i);
    end
    MtM_mu2I_inv = 1./(MtM + mu2 + eps);
    mu1_EtE_mu2_DtHtHD_inv = 1./(mu1*EtE + mu2*Dt.*sum(abs(H).^2,4).*D + eps);

    
    %% Allocate memory and initialize
    
    % Dual and aux variables for Gx = z1
    z1_1 = zeros(N); z1_2 = zeros(N); z1_3 = zeros(N);
    s1_1 = zeros(N); s1_2 = zeros(N); s1_3 = zeros(N);
    
    % Dual and aux variables for HDFx = z2
    s2 = zeros([N,numKernels]); z2 = zeros([N,numKernels]);
    
    Fx_prev = zeros(N);
     
    %% Iteratively solve
    tic
    for outerIter = 1:maxOuterIter
        
        % Update x
        Fz1_minus_s1_1 = fftn(z1_1 - s1_1); % x direction
        Fz1_minus_s1_2 = fftn(z1_2 - s1_2); % y direction
        Fz1_minus_s1_3 = fftn(z1_3 - s1_3); % z direction
        
        Fx = mu1_EtE_mu2_DtHtHD_inv.* (mu1*(E1t.*Fz1_minus_s1_1 + E2t.*Fz1_minus_s1_2 + E3t.*Fz1_minus_s1_3) + mu2*Dt.*sum(Ht.*(z2-s2),4));
        
        % Update aux variable z2: z2 = HDFx
        HDFx = bsxfun(@times,H,D.*Fx);
        for i = 1:numKernels
            z2(:,:,:,i) = fftn(MtM_mu2I_inv(:,:,:,i).*(MtM_FtHF_phi(:,:,:,i) + mu2*ifftn(HDFx(:,:,:,i) + s2(:,:,:,i))));
        end

        % Update aux variable z1: z1 = Gx
        Gx_1 = ifftn(E1.*Fx); Gx_2 = ifftn(E2.*Fx); Gx_3 = ifftn(E3.*Fx); % Compute gradient along x, y, and z
    	
        z1_1 = max(abs(Gx_1 + s1_1)-alpha/mu1,0).*sign(Gx_1 + s1_1);
        z1_2 = max(abs(Gx_2 + s1_2)-alpha/mu1,0).*sign(Gx_2 + s1_2);
        z1_3 = max(abs(Gx_3 + s1_3)-alpha/mu1,0).*sign(Gx_3 + s1_3);

        % Update the scaled dual variables s1, s2
        s1_1 = s1_1 + Gx_1 - z1_1;
        s1_2 = s1_2 + Gx_2 - z1_2;
        s1_3 = s1_3 + Gx_3 - z1_3;
        
        s2 = s2 + HDFx - z2;
        
        % Stopping criterion
        if outerIter > 1
            rmse_soln = 100*norm(Fx_prev(:)-Fx(:))/(eps+norm(Fx_prev(:)));
            %disp(['Iter: ', num2str(outerIter), '  Update: ', num2str(rmse_soln)])
            
            if rmse_soln < tol_soln 
                break
            end
        end
        Fx_prev = Fx;
         
    end
    toc
    
    out.x = real(ifftn(Fx));
    out.rmse_soln = rmse_soln;
end
