clear, clc
set(0,'DefaultFigureWindowStyle','docked')
data_path = [pwd,'/Data/numerical_qsm_phantom'];
load([data_path,'/numerical_qsm_phantom.mat'])
stop
N = size(chi_total);
TE = 21e-3;             % second
B0 = 2.89;              % Tesla
gyro = 2*pi*42.58;

% For display
display_position = round(cropsize/2);
display_orientation = [90,90,90];
display_caxis = [-0.1,0.1];
display_caxis_diff = display_caxis/5;

%% dipole kernel
D = create_dipole_kernel(B0_dir, voxel_size, N, 1);

%% total and tissue phase
load([data_path,'/noise_2pt5_percent'])
stop
phase_total_clean = real(ifftn(D .* fftn(chi_total))); 
phase_total = real(ifftn(D .* fftn(chi_total))) + noise; 
rmse_from_noise = 100*norm(wrapToPi(phase_total_clean(:))-wrapToPi(phase_total(:)))/norm(wrapToPi(phase_total_clean(:)));
disp(['RMSE from noise = ', num2str(rmse_from_noise)])
phase_tissue = real(ifftn(D .* fftn(chi_tissue))); 

imagesc3d2(phase_total(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1), display_position, 2, display_orientation, [-.1,.1], [], 'Total phase')

%% Generate SMV kernels and masks
min_radius = 1;
max_radius = 5;
step_size_radius = 1;

out = create_SMVkernel(phase_total, chi_mask, min_radius, max_radius, step_size_radius, N, voxel_size);

SMV_kernels = out.SMV_kernel;
SMV_inv_kernel = out.SMV_inv_kernel;
SMV_masks = out.SMV_mask;
SMV_phase = out.SMV_phase;
mask_Sharp = out.mask_eval;

%% Ground truth
chi_tissue_0mean = zeros(N);
chi_tissue_0mean(mask_Sharp==1) = chi_tissue(mask_Sharp==1) - mean(chi_tissue(mask_Sharp==1));
chi_tissue_0mean_disp = chi_tissue_0mean - (mask_Sharp == 0);
imagesc3d2(chi_tissue_0mean_disp(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1), display_position, 1, display_orientation, display_caxis, [], 'Underlying susceptibility map')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Single-Step QSM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Single-step L2 QSM
params = [];
params.alpha = 3e-4;                % Regularization param for quadratic smoothing
params.N = N;                       % Number of voxels
params.M = SMV_masks;               % Mask for each reliable region
params.H = SMV_kernels;             % SMV kernels
params.D = D;                       % Dipole kernel
params.phase_unwrap = phase_total;  % Unwrapped total phase

out_ss_l2 = SS_L2_QSM(params);

chi_ss_l2 = mask_Sharp .* out_ss_l2.x;

chi_ss_l2_0mean = zeros(N);
chi_ss_l2_0mean(mask_Sharp==1) = chi_ss_l2(mask_Sharp==1) - mean(chi_ss_l2(mask_Sharp==1));

rmse_ss_l2 = 100 * norm((chi_ss_l2_0mean(:) - chi_tissue_0mean(:)) .* mask_Sharp(:)) / norm(chi_tissue_0mean(:) .* mask_Sharp(:));
disp(['RMSE Single-Step L2 = ', num2str(rmse_ss_l2)])

% Display results
chi_disp = chi_ss_l2_0mean - (mask_Sharp==0);
imagesc3d2(chi_disp(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1), display_position, 11, display_orientation, display_caxis, [], ['Single-Step L2: ', num2str(rmse_ss_l2)])

chi_diff_disp = chi_ss_l2_0mean - chi_tissue_0mean - (mask_Sharp == 0);
imagesc3d2(chi_diff_disp(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1), display_position, 21, display_orientation, display_caxis_diff, [], ['Single-Step L2: ',num2str(rmse_ss_l2)])

%% Single-step TV QSM
params = [];
params.alpha = 7e-5;                % Regularization param for ||Gx||_1
params.mu1 = 3e-2;                  % Aug Lagrangian param for z1 = Gx
params.mu2 = params.mu1;            % Aug Lagrangian param for z2 = HDFx
params.maxOuterIter = 100;          % Max number of iter
params.tol_soln = 1;                % Stopping criterion: RMSE change in solution 
params.N = N;                       % Number of voxels
params.M = SMV_masks;               % RMask for each reliable region
params.H = SMV_kernels;             % SMV kernels
params.D = D;                       % Dipole kernel
params.phase_unwrap = phase_total;  % Unwrapped total phase
 
out_ss_tv = SS_TV_QSM(params);

chi_ss_tv = mask_Sharp .* out_ss_tv.x;

chi_ss_tv_0mean = zeros(N);
chi_ss_tv_0mean(mask_Sharp==1) = chi_ss_tv(mask_Sharp==1) - mean(chi_ss_tv(mask_Sharp==1));

rmse_ss_tv = 100 * norm((chi_ss_tv_0mean(:) - chi_tissue_0mean(:)) .* mask_Sharp(:)) / norm(chi_tissue_0mean(:) .* mask_Sharp(:));
disp(['RMSE Single-Step TV = ', num2str(rmse_ss_tv)])

% Display results
chi_disp = chi_ss_tv_0mean  - (mask_Sharp==0);
imagesc3d2(chi_disp(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1), display_position, 12, display_orientation, display_caxis, [], ['Single-Step TV: ',num2str(rmse_ss_tv)])

chi_diff_disp = chi_ss_tv_0mean - chi_tissue_0mean - (mask_Sharp == 0);
imagesc3d2(chi_diff_disp(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1), display_position, 22, display_orientation, display_caxis_diff, [], ['Single-Step TV: ',num2str(rmse_ss_tv)])

%% Single-step TGV QSM
params = [];
alpha1 = 6e-5;
params.alpha0 = 2*alpha1;           % Regularization param for ||Ev||_1
params.alpha1 = alpha1;             % Regularization param for ||Gx-v||_1
params.mu0 = 3e-2;                  % Aug Lagrangian param for z0 = Ev
params.mu1 = params.mu0;            % Aug Lagrangian param for z1 = Gx
params.mu2 = params.mu0;            % Aug Lagrangian param for z2 = HDFx
params.maxOuterIter = 100;          % Max number of iter
params.tol_soln = 1;                % Stopping criterion: RMSE change in solution 
params.N = N;                       % Number of voxels
params.M = SMV_masks;               % Mask for each reliable region
params.H = SMV_kernels;             % SMV kernels
params.D = D;                       % Dipole kernel
params.phase_unwrap = phase_total;  % Unwrapped total phase

out_ss_tgv = SS_TGV_QSM(params);

chi_ss_tgv = mask_Sharp .* out_ss_tgv.x;

chi_ss_tgv_0mean = zeros(N);
chi_ss_tgv_0mean(mask_Sharp==1) = chi_ss_tgv(mask_Sharp==1) - mean(chi_ss_tgv(mask_Sharp==1));

rmse_ss_tgv = 100 * norm((chi_ss_tgv_0mean(:) - chi_tissue_0mean(:)) .* mask_Sharp(:)) / norm(chi_tissue_0mean(:) .* mask_Sharp(:));
disp(['RMSE Single-Step TGV = ', num2str(rmse_ss_tgv)])

% Display results
chi_disp = chi_ss_tgv_0mean  - (mask_Sharp==0);
imagesc3d2(chi_disp(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1), display_position, 13, display_orientation, display_caxis, [], ['Single-Step TGV: ',num2str(rmse_ss_tgv)])

chi_diff_disp = chi_ss_tgv_0mean - chi_tissue_0mean - (mask_Sharp == 0);
imagesc3d2(chi_diff_disp(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1), display_position, 23, display_orientation, display_caxis_diff, [], ['Single-Step TGV: ',num2str(rmse_ss_tgv)])

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Multiple-Step QSM %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inverse V-SHARP
SMV_thres = 2e-2;

SMV_inv_kernel_temp = zeros(N);
SMV_inv_kernel_temp( abs(SMV_inv_kernel) > SMV_thres ) = 1 ./ SMV_inv_kernel( abs(SMV_inv_kernel) > SMV_thres );

phase_inv_Sharp = mask_Sharp .* ifftn(SMV_inv_kernel_temp .* fftn(SMV_phase));

rmse_inv_Sharp = 100 * norm((phase_inv_Sharp(:) - phase_tissue(:)) .* mask_Sharp(:)) / norm(phase_tissue(:) .* mask_Sharp(:));
disp(['RMSE inverse Sharp = ', num2str(rmse_inv_Sharp)])

%% Multiple-step L2 QSM
params = [];
params.alpha = 2e-2;                    % Regularization param for quadratic smoothing
params.N = N;                           % Number of voxels
params.D = D;                           % Dipole kernel
params.phase_unwrap = phase_inv_Sharp;  % V-SHARP filtered phase

out_ms_l2 = MS_L2_QSM(params);

chi_ms_l2 = out_ms_l2.x .* mask_Sharp;
chi_ms_l2_0mean = zeros(N);
chi_ms_l2_0mean(mask_Sharp==1) = chi_ms_l2(mask_Sharp==1) - mean(chi_ms_l2(mask_Sharp==1));

rmse_ms_l2 = 100 * norm((chi_ms_l2_0mean(:) - chi_tissue_0mean(:)) .* mask_Sharp(:)) / norm(chi_tissue_0mean(:) .* mask_Sharp(:));
disp(['RMSE VSHARP L2 = ', num2str(rmse_ms_l2)])

% Display results
chi_disp = chi_ms_l2_0mean  - (mask_Sharp==0);
imagesc3d2(chi_disp(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1), display_position, 14, display_orientation, display_caxis, [], ['VSHARP L2: ', num2str(rmse_ms_l2)])

chi_diff_disp = chi_ms_l2_0mean - chi_tissue_0mean - (mask_Sharp == 0);
imagesc3d2(chi_diff_disp(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1), display_position, 24, display_orientation, display_caxis_diff, [], ['VSHARP L2: ',num2str(rmse_ms_l2)])

%% Multiple-step TV QSM
params = [];
params.alpha = 2e-4;                    % Regularization param for ||Gx||_1
params.mu = 3e-2;                       % Aug Lagrangrian param for z = Gx
params.maxOuterIter = 50;               % Max number of iterations
params.tol_soln = 1;                    % Stopping criterion: RMSE change in solution 
params.N = N;                           % Number of voxels
params.D = D;                           % Dipole kernel
params.phase_unwrap = phase_inv_Sharp;  % V-SHARP filtered phase

out_ms_tv = MS_TV_QSM(params);

chi_ms_tv = out_ms_tv.x .* mask_Sharp;
chi_ms_tv_0mean = zeros(N);
chi_ms_tv_0mean(mask_Sharp==1) = chi_ms_tv(mask_Sharp==1) - mean(chi_ms_tv(mask_Sharp==1));

rmse_ms_tv = 100 * norm((chi_ms_tv_0mean(:) - chi_tissue_0mean(:)) .* mask_Sharp(:)) / norm(chi_tissue_0mean(:) .* mask_Sharp(:));
disp(['RMSE VSHARP TV = ', num2str(rmse_ms_tv)])

% Display results
chi_disp = chi_ms_tv_0mean  - (mask_Sharp==0);
imagesc3d2(chi_disp(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1), display_position, 15, display_orientation, display_caxis, [], ['VSHARP TV: ', num2str(rmse_ms_tv)])

chi_diff_disp = chi_ms_tv_0mean - chi_tissue_0mean - (mask_Sharp == 0);
imagesc3d2(chi_diff_disp(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1), display_position, 25, display_orientation, display_caxis_diff, [], ['VSHARP TV: ',num2str(rmse_ms_tv)])

%% Multiple-step TGV QSM
params = [];
params.alpha1 = 2e-4;                   % Regularization param for ||Gx-v||_1
params.alpha0 = 2 * params.alpha1;      % Regularization param for ||Ev||_1
params.mu1 = 3e-2;                      % Aug Lagrangrian param for z1 = Gx-v
params.mu0 = params.mu1;                % Aug Lagrangrian param for z0 = Ev
params.maxOuterIter = 50;               % Max number of iterations
params.tol_update = 1;                  % Stopping criterion: RMSE change in solution 
params.N = N;                           % Number of voxels
params.kspace = fftn(phase_inv_Sharp);  % DTFT of V-SHARP filtered phase
params.D = D;                           % Dipole kernel
 
out_ms_tgv = MS_TGV_QSM(params); 

chi_ms_tgv = out_ms_tgv.x .* mask_Sharp;
chi_ms_tgv_0mean = zeros(N);
chi_ms_tgv_0mean(mask_Sharp==1) = chi_ms_tgv(mask_Sharp==1) - mean(chi_ms_tgv(mask_Sharp==1));

rmse_ms_tgv = 100 * norm((chi_ms_tgv_0mean(:) - chi_tissue_0mean(:)) .* mask_Sharp(:)) / norm(chi_tissue_0mean(:) .* mask_Sharp(:));
disp(['RMSE VSHARP TGV = ', num2str(rmse_ms_tgv)])

% Display results
chi_disp = chi_ms_tgv_0mean  - (mask_Sharp==0);
imagesc3d2(chi_disp(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1), display_position, 16, display_orientation, display_caxis, [], ['VSHARP TGV: ', num2str(rmse_ms_tgv)])

chi_diff_disp = chi_ms_tgv_0mean - chi_tissue_0mean - (mask_Sharp == 0);
imagesc3d2(chi_diff_disp(1:d1last-d1first+1, 1:d2last-d2first+1, 1:d3last-d3first+1), display_position, 26, display_orientation, display_caxis_diff, [], ['VSHARP TGV: ',num2str(rmse_ms_tgv)])

