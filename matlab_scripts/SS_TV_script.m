function [chi_ss_tv_ppb]=SS_TV_script(freqLoc,maskLoc,reliabilityMaskLoc,CF,alpha,B0_dir,SuscFilename)    
    addpath([dirname(mfilename('fullpath')),'/TGV_SS_QSM/SS_TGV_QSM_Toolbox']);
    addpath([dirname(mfilename('fullpath')),'/NIfTI_20140122']);
   
    CF=str2double(CF);
    alpha=str2double(alpha);
    B0_dir=str2double(B0_dir);
  
    %mask
    chi_mask_obj=load_untouch_nii(maskLoc);
    chi_mask=double(chi_mask_obj.img);
    
    %reliability
    reliability_mask_obj=load_untouch_nii(reliabilityMaskLoc);
    reliability_mask=double(reliability_mask_obj.img);

    %fieldmap
    chi_total_obj=load_untouch_nii(freqLoc);
    phase_total=double(chi_total_obj.img);
        
    if B0_dir==1
        %B0_dir=[1;0;0];
        %put B0 in 3rd dimension
        chi_mask=permute(chi_mask,[3,2,1]);
        reliability_mask=permute(reliability_mask,[3,2,1]);
        phase_total=permute(phase_total,[3,2,1]);
    elseif B0_dir==2
        %B0_dir=[0;1;0];
        %put B0 in 3rd dimension
        chi_mask=permute(chi_mask,[1,3,2]);
        reliability_mask=permute(reliability_mask,[1,3,2]);
        phase_total=permute(phase_total,[1,3,2]);
    elseif B0_dir==3
        %B0_dir=[0;0;1];
    end
    
    aktmp=size(phase_total);  
    
    if mod(aktmp(end),2)==0
        phase_total=padarray(phase_total,[0,0,10],0,'pre');
        chi_mask=padarray(chi_mask,[0,0,10],0,'pre');
        reliability_mask=padarray(reliability_mask,[0,0,10],0,'pre');
    else
        phase_total=padarray(phase_total,[0,0,11],0,'pre');
        chi_mask=padarray(chi_mask,[0,0,11],0,'pre');
        reliability_mask=padarray(reliability_mask,[0,0,11],0,'pre');
    end
    
    %voxel_size
    voxel_size=chi_total_obj.hdr.dime.pixdim(2:4);
    
    N = size(phase_total);

    %% dipole kernel
    D = create_dipole_kernel(B0_dir, voxel_size, N, 1);

    %% Generate SMV kernels and masks
    min_radius = 1;
    max_radius = 5;
    step_size_radius = 1;

    out = create_SMVkernel(phase_total, chi_mask, min_radius, max_radius, step_size_radius, N, voxel_size);

    SMV_kernels = out.SMV_kernel;
    SMV_masks = out.SMV_mask.*reliability_mask; 
    mask_Sharp = out.mask_eval;   

    %% Single-step TV QSM
    params = [];
    params.alpha = alpha;               % Regularization param for ||Gx||_1
    params.mu1 = 0.1;%3e-2;                  % Aug Lagrangian param for z1 = Gx
    params.mu2 = params.mu1;            % Aug Lagrangian param for z2 = HDFx
    %!! note: lower number of iterations in SB doesn't seem to affect Lcurve
    params.maxOuterIter = 200;          % Max number of iter
    %!! note: lower number of iterations in SB doesn't seem to affect Lcurve
    params.tol_soln = 0.1;                % Stopping criterion: RMSE change in solution 
    params.N = N;                       % Number of voxels
    params.M = SMV_masks;               % RMask for each reliable region
    params.H = SMV_kernels;             % SMV kernels
    params.D = D;                       % Dipole kernel
    params.phase_unwrap = phase_total;  % Unwrapped total phase




out_ss_tv = SS_TV_QSM(params);
chi_ss_tv = mask_Sharp .* out_ss_tv.x;

chi_ss_tv_0mean = zeros(N);
chi_ss_tv_0mean(mask_Sharp==1) = chi_ss_tv(mask_Sharp==1) - mean(chi_ss_tv(mask_Sharp==1));
chi_ss_tv_ppb=chi_ss_tv_0mean/((CF*1e-9)*2*pi);

if mod(aktmp(end),2)==0
    chi_ss_tv_ppb=chi_ss_tv_ppb(:,:,11:end); 
else
    chi_ss_tv_ppb=chi_ss_tv_ppb(:,:,12:end);
end

if B0_dir==1
    %B0_dir=[1;0;0];
    %put B0 back in 1st dimension
    chi_ss_tv_ppb=permute(chi_ss_tv_ppb,[3,2,1]);      
elseif B0_dir==2
    %B0_dir=[0;1;0];
    %put B0 back in 2nd dimension
    chi_ss_tv_ppb=permute(chi_ss_tv_ppb,[1,3,2]);
elseif B0_dir==3
    %B0_dir=[0;0;1];
end

tmp=make_nii(chi_ss_tv_ppb);
tmp.hdr=chi_total_obj.hdr;
save_nii(tmp,SuscFilename);
return



