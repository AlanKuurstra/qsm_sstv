function out = MS_TV_QSM(params)

% Retrieve data
mu = params.mu;             % gradient consistency
alpha = params.alpha;     % gradient L1 penalty
maxOuterIter = params.maxOuterIter;
tol_soln = params.tol_soln;
N = params.N;
D = params.D;
phase_unwrap = params.phase_unwrap;

z_dx = zeros(N);      z_dy = zeros(N);      z_dz = zeros(N);
s_dx = zeros(N);      s_dy = zeros(N);      s_dz = zeros(N);

x = zeros(N);

kspace = fftn(phase_unwrap);
Dt_kspace = conj(D) .* kspace;

[kx, ky, kz] = ndgrid(0:N(1)-1, 0:N(2)-1, 0:N(3)-1);
Ex = 1 - exp(2i .* pi .* kx / N(1));
Ey = 1 - exp(2i .* pi .* ky / N(2));
Ez = 1 - exp(2i .* pi .* kz / N(3));

Ext = conj(Ex);     Eyt = conj(Ey);     Ezt = conj(Ez);
EtE = Ext .* Ex + Eyt .* Ey + Ezt .* Ez;
DtD = abs(D).^2;

tic
for t = 1:maxOuterIter
    % update x : susceptibility estimate
    tx = Ext .* fftn(z_dx - s_dx);
    ty = Eyt .* fftn(z_dy - s_dy);
    tz = Ezt .* fftn(z_dz - s_dz);
    
    x_prev = x;
    x = ifftn( (mu * (tx + ty + tz) + Dt_kspace) ./ (eps + DtD + mu * EtE) );

    x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
%     disp(['Iter: ', num2str(t), '   Update: ', num2str(x_update)])
    
    if x_update < tol_soln
        break
    end
    
    if t < maxOuterIter
        % update z : gradient varible
        Fx = fftn(x);
        x_dx = ifftn(Ex .* Fx);
        x_dy = ifftn(Ey .* Fx);
        x_dz = ifftn(Ez .* Fx);

        z_dx = max(abs(x_dx + s_dx) - alpha / mu, 0) .* sign(x_dx + s_dx);
        z_dy = max(abs(x_dy + s_dy) - alpha / mu, 0) .* sign(x_dy + s_dy);
        z_dz = max(abs(x_dz + s_dz) - alpha / mu, 0) .* sign(x_dz + s_dz);

        % update s : Lagrange multiplier
        s_dx = s_dx + x_dx - z_dx;
        s_dy = s_dy + x_dy - z_dy;            
        s_dz = s_dz + x_dz - z_dz;            
    end
end
toc
out.x = real(x);

end

