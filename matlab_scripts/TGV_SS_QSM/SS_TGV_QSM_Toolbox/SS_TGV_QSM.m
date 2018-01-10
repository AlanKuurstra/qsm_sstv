function out = SS_TGV_QSM(params)

    % Retrieve data
    alpha0 = params.alpha0;
    alpha1 = params.alpha1;
    mu0 = params.mu0;
    mu1 = params.mu1;
    mu2 = params.mu2;
    maxOuterIter = params.maxOuterIter;
    tol_soln = params.tol_soln;
    N = params.N;
    M = params.M;
    H = params.H;
    D = params.D;
    numKernels = size(H,4);

    % Precompute
    
    [k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);
    E1 = 1 - exp(2i .* pi .* k1 / N(1)); % G1 = F'*E1*F
    E2 = 1 - exp(2i .* pi .* k2 / N(2)); % G2 = F'*E2*F
    E3 = 1 - exp(2i .* pi .* k3 / N(3)); % G3 = F'*E3*F

    E1t = conj(E1); E2t = conj(E2); E3t = conj(E3);
    
    E1tE1 = E1t.*E1; E2tE2 = E2t.*E2; E3tE3 = E3t.*E3;
    mu0_over_2_E1tE2 = mu0/2*E1t.*E2;
    mu0_over_2_E1tE3 = mu0/2*E1t.*E3;
    mu0_over_2_E2tE3 = mu0/2*E2t.*E3;
    
    Dt = conj(D); Ht = conj(H);

    mu1_E_sos_mu2_DtHtHD    = mu1*(E1tE1 + E2tE2 + E3tE3) + mu2*Dt.*sum(abs(H).^2,4).*D;
    mu1I_mu0_E_wsos1        = mu1 + mu0*(E1tE1 + (E2tE2 + E3tE3)/2);
    mu1I_mu0_E_wsos2        = mu1 + mu0*(E1tE1/2 + E2tE2 + E3tE3/2);
    mu1I_mu0_E_wsos3        = mu1 + mu0*((E1tE1 + E2tE2)/2 + E3tE3);
    
    MtM = abs(M).^2;
    MtM_FtHF_phi = zeros([N,numKernels]); F_phi = fftn(params.phase_unwrap);
    for i = 1:numKernels
        MtM_FtHF_phi(:,:,:,i) = MtM(:,:,:,i).*ifftn(H(:,:,:,i).*F_phi);
    end
    MtM_mu2I_inv = 1./(MtM + mu2 + eps);
    
    % Precompute for Cramer's rule (we can also avoid repetitively computing things in here...) 
    
    % A = [a1, a5t, a6t, a8t; a5, a2, a7t, a9t; a6, a7, a3, a10t; a8, a9,a10, a4]
    a1 = mu1_E_sos_mu2_DtHtHD; a2 = mu1I_mu0_E_wsos1; a3 = mu1I_mu0_E_wsos2; a4 = mu1I_mu0_E_wsos3;
    a5 = -mu1*E1; a6 = -mu1*E2; a7 = mu0_over_2_E1tE2; a8 = -mu1*E3; a9 = mu0_over_2_E1tE3; a10 = mu0_over_2_E2tE3;
    a5t = conj(a5); a6t = conj(a6); a7t = conj(a7); a8t = conj(a8); a9t = conj(a9); a10t = conj(a10);    
    
    % For x
    D11 = a2.*a3.*a4 + a7t.*a9.*a10t + a7.*a9t.*a10 - a3.*a9.*a9t - a2.*a10.*a10t - a4.*a7.*a7t;
    D21 = a3.*a4.*a5t  + a6t.*a9.*a10t + a7.*a8t.*a10  - a3.*a8t.*a9 - a5t.*a10.*a10t - a4.*a6t.*a7;
    D31 = a4.*a5t.*a7t + a6t.*a9.*a9t  + a2.*a8t.*a10  - a7t.*a8t.*a9 - a5t.*a9t.*a10 - a2.*a4.*a6t;
    D41 = a5t.*a7t.*a10t + a6t.*a7.*a9t + a2.*a3.*a8t     - a7.*a7t.*a8t - a3.*a5t.*a9t    - a2.*a6t.*a10t;

    % For vx
    D12 = a3.*a4.*a5 + a7t.*a8.*a10t + a6.*a9t.*a10 - a3.*a8.*a9t - a5.*a10.*a10t - a4.*a6.*a7t;
    D22 = a1.*a3.*a4 + a6t.*a8.*a10t + a6.*a8t.*a10 - a3.*a8.*a8t - a1.*a10.*a10t - a4.*a6.*a6t;
    D32 = a1.*a4.*a7t + a6t.*a8.*a9t + a5.*a8t.*a10 - a7t.*a8.*a8t - a1.*a9t.*a10 - a4.*a5.*a6t;
    D42 = a1.*a7t.*a10t + a6.*a6t.*a9t + a3.*a5.*a8t - a6.*a7t.*a8t - a1.*a3.*a9t - a5.*a6t.*a10t;

    % For vy
    D13 = a4.*a5.*a7 + a2.*a8.*a10t + a6.*a9.*a9t - a7.*a8.*a9t - a5.*a9.*a10t - a2.*a4.*a6;
    D23 = a1.*a4.*a7 + a5t.*a8.*a10t +a6.*a8t.*a9 - a7.*a8.*a8t - a1.*a9.*a10t - a4.*a5t.*a6;
    D33 = a1.*a2.*a4 + a5t.*a8.*a9t + a5.*a8t.*a9 - a2.*a8.*a8t - a1.*a9.*a9t - a4.*a5.*a5t;
    D43 = a1.*a2.*a10t + a5t.*a6.*a9t + a5.*a7.*a8t - a2.*a6.*a8t - a1.*a7.*a9t - a5.*a5t.*a10t;

    % For vz
    D14 = a5.*a7.*a10 + a2.*a3.*a8 + a6.*a7t.*a9 - a7.*a7t.*a8 - a3.*a5.*a9 -a2.*a6.*a10;
    D24 = a1.*a7.*a10 + a3.*a5t.*a8 + a6.*a6t.*a9 - a6t.*a7.*a8 - a1.*a3.*a9 - a5t.*a6.*a10;
    D34 = a1.*a2.*a10 + a5t.*a7t.*a8 + a5.*a6t.*a9 - a2.*a6t.*a8 - a1.*a7t.*a9 - a5.*a5t.*a10;
    D44 = a1.*a2.*a3 + a5t.*a6.*a7t + a5.*a6t.*a7 - a2.*a6.*a6t - a1.*a7.*a7t - a3.*a5.*a5t;

    det_A =     (a1.*D11 - a5.*D21 + a6.*D31 - a8.*D41) + eps;
    
    % Allocate memory and initialize
    
    % Dual and aux variables for Ev = z0
    s0_1 = zeros(N); z0_1 = zeros(N); % xx
    s0_2 = zeros(N); z0_2 = zeros(N); % yy
    s0_3 = zeros(N); z0_3 = zeros(N); % zz
    s0_4 = zeros(N); z0_4 = zeros(N); % xy
    s0_5 = zeros(N); z0_5 = zeros(N); % xz
    s0_6 = zeros(N); z0_6 = zeros(N); % yz
    
    % Dual and aux variables for Gx = z1
    s1_1 = zeros(N); z1_1 = zeros(N); % x
    s1_2 = zeros(N); z1_2 = zeros(N); % y
    s1_3 = zeros(N); z1_3 = zeros(N); % z
    
    % Dual and aux variables for HDFx = z2
    s2 = zeros([N,numKernels]); z2 = zeros([N,numKernels]);
    
    Fx_prev = zeros(N);
     
    % Iteratively solve
    tic
    for outerIter = 1:maxOuterIter
        
        % Update x
        F_z0_minus_s0_1 = fftn(z0_1 - s0_1); % xx
        F_z0_minus_s0_2 = fftn(z0_2 - s0_2); % yy
        F_z0_minus_s0_3 = fftn(z0_3 - s0_3); % zz
        F_z0_minus_s0_4 = fftn(z0_4 - s0_4); % xy
        F_z0_minus_s0_5 = fftn(z0_5 - s0_5); % xz
        F_z0_minus_s0_6 = fftn(z0_6 - s0_6); % yz

        F_z1_minus_s1_1 = fftn(z1_1 - s1_1); % x
        F_z1_minus_s1_2 = fftn(z1_2 - s1_2); % y
        F_z1_minus_s1_3 = fftn(z1_3 - s1_3); % z
        
        % Compute the right hand side
        rhs1    =  mu2*Dt.*sum(Ht.*(z2-s2),4) + mu1*(E1t.*F_z1_minus_s1_1 + E2t.*F_z1_minus_s1_2 + E3t.*F_z1_minus_s1_3);
        rhs2    = -mu1*F_z1_minus_s1_1       + mu0*(E1t.*F_z0_minus_s0_1 + E2t.*F_z0_minus_s0_4 + E3t.*F_z0_minus_s0_5);
        rhs3    = -mu1*F_z1_minus_s1_2       + mu0*(E2t.*F_z0_minus_s0_2 + E1t.*F_z0_minus_s0_4 + E3t.*F_z0_minus_s0_6);
        rhs4    = -mu1*F_z1_minus_s1_3       + mu0*(E3t.*F_z0_minus_s0_3 + E1t.*F_z0_minus_s0_5 + E2t.*F_z0_minus_s0_6);
        
        % Simultaneously apply Cramer's rule to the 4x4 matrices
        Fx = (rhs1.*D11 -rhs2.*D21 +rhs3.*D31 - rhs4.*D41)./det_A;
        Fv1 = (-rhs1.*D12 +rhs2.*D22 -rhs3.*D32 + rhs4.*D42)./det_A;
        Fv2 = (rhs1.*D13 -rhs2.*D23 +rhs3.*D33 - rhs4.*D43)./det_A;
        Fv3 = (-rhs1.*D14 +rhs2.*D24 -rhs3.*D34 + rhs4.*D44)./det_A;  
        
        v1 = ifftn(Fv1);
        v2 = ifftn(Fv2);
        v3 = ifftn(Fv3);
        
        % Update aux variable z2: z2 = HDFx
        HDFx = bsxfun(@times,H,D.*Fx);
        for i = 1:numKernels
            z2(:,:,:,i) = fftn(MtM_mu2I_inv(:,:,:,i).*(MtM_FtHF_phi(:,:,:,i) + mu2*ifftn(HDFx(:,:,:,i) + s2(:,:,:,i))));
        end
        
        % Compute gradients for z0 and z1 update
        Ev_1 = ifftn(E1.*Fv1);
        Ev_2 = ifftn(E2.*Fv2);
        Ev_3 = ifftn(E3.*Fv3);
        Ev_4 = ifftn(E1.*Fv2 + E2.*Fv1)/2;
        Ev_5 = ifftn(E1.*Fv3 + E3.*Fv1)/2;
        Ev_6 = ifftn(E2.*Fv3 + E3.*Fv2)/2;
        
        Gx_1 = ifftn(E1.*Fx);
        Gx_2 = ifftn(E2.*Fx);
        Gx_3 = ifftn(E3.*Fx);
        
        % Update aux variable z0: z0 = Ev
        z0_1 = max(abs(Ev_1 + s0_1)-alpha0/mu0,0).*sign(Ev_1 + s0_1);
        z0_2 = max(abs(Ev_2 + s0_2)-alpha0/mu0,0).*sign(Ev_2 + s0_2);
        z0_3 = max(abs(Ev_3 + s0_3)-alpha0/mu0,0).*sign(Ev_3 + s0_3);
        z0_4 = max(abs(Ev_4 + s0_4)-alpha0/mu0,0).*sign(Ev_4 + s0_4);
        z0_5 = max(abs(Ev_5 + s0_5)-alpha0/mu0,0).*sign(Ev_5 + s0_5);
        z0_6 = max(abs(Ev_6 + s0_6)-alpha0/mu0,0).*sign(Ev_6 + s0_6);
        
        % Update aux variable z1: z1 = Gx
        Gx_minus_v_1 = Gx_1 - v1; Gx_minus_v_2 = Gx_2 - v2; Gx_minus_v_3 = Gx_3 - v3;
        z1_1 = max(abs(Gx_minus_v_1 + s1_1)-alpha1/mu1,0).*sign(Gx_minus_v_1 + s1_1);
        z1_2 = max(abs(Gx_minus_v_2 + s1_2)-alpha1/mu1,0).*sign(Gx_minus_v_2 + s1_2);
        z1_3 = max(abs(Gx_minus_v_3 + s1_3)-alpha1/mu1,0).*sign(Gx_minus_v_3 + s1_3);

        % Update the scaled dual variables s0, s1, and  s2
        s0_1 = s0_1 + Ev_1 - z0_1;
        s0_2 = s0_2 + Ev_2 - z0_2;
        s0_3 = s0_3 + Ev_3 - z0_3;
        s0_4 = s0_4 + Ev_4 - z0_4;
        s0_5 = s0_5 + Ev_5 - z0_5;
        s0_6 = s0_6 + Ev_6 - z0_6;
        
        s1_1 = s1_1 + Gx_minus_v_1 - z1_1;
        s1_2 = s1_2 + Gx_minus_v_2 - z1_2;
        s1_3 = s1_3 + Gx_minus_v_3 - z1_3;
        
        s2 = s2 + HDFx - z2;
        
        % Stopping criterion
        if outerIter > 1
            rmse_soln = 100*norm(Fx_prev(:)-Fx(:))/(eps+norm(Fx_prev(:)));
%             disp(['Iter: ', num2str(outerIter), '  Update: ', num2str(rmse_soln)])
            
            if rmse_soln < tol_soln 
                break
            end
        end
        Fx_prev = Fx;
        
    end
    toc
    
    out.x = real(ifftn(Fx));
    out.v1 = v1;
    out.v2 = v2;
    out.v3 = v3;
    out.rmse_soln = rmse_soln;
end

